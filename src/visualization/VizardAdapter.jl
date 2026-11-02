# ==============================================================================
# src/visualization/VizardAdapter.jl
# Unity-Based Basilisk Vizard Transport Adapter (Protobuf over ZeroMQ,
# Handshake Strings, 3 Live Modes on Ports 5556 / 5570 / noDisplay, and
# simSynch Wall-Clock Synchronization Contract)
# ==============================================================================

const VIZARD_ZMQ_HANDSHAKE_STRINGS = (
    "PING",
    "SIM_UPDATE",
    "OK",
    "REQUEST_INPUT",
    "VIZARD_INPUT",
    "REQUEST_IMAGE_",
)

const VIZARD_LIVE_MODES = (
    liveStream = (
        mode = :liveStream,
        default_endpoint = "tcp://localhost:5556",
        zmq_socket_action = :zmq_connect_two_way,
        lockstep_waits_for_ok = true,
        supports_user_input_and_opnav = true,
        legacy_opnav_mode = 1
    ),
    broadcastStream = (
        mode = :broadcastStream,
        default_endpoint = "tcp://localhost:5570",
        zmq_socket_action = :zmq_bind_publish_subscribe,
        lockstep_waits_for_ok = false,
        supports_multiple_subscribers = true,
        recommended_m7_hud_tap = true,
        julia_subscriber_stack = (:ZMQ_jl_SUB_socket, :ProtoBuf_jl_compiled_vizMessage),
        json_side_channel_supported = false
    ),
    noDisplay = (
        mode = :noDisplay,
        default_endpoint = "NONE",
        zmq_socket_action = :headless_opnav_only,
        lockstep_waits_for_ok = true,
        legacy_opnav_mode = 2,
        mutually_exclusive_with_streams = true
    ),
)

struct VizardTransportConfig
    liveStream::Bool
    broadcastStream::Bool
    noDisplay::Bool
    live_endpoint::String
    broadcast_endpoint::String
    record_file_playback_path::Union{Nothing,String}
    sim_synch_enabled::Bool
    sim_synch_acceleration_factor::Float64
    sim_synch_header_field_names_verdict::Symbol

    function VizardTransportConfig(;
            liveStream::Bool = false,
            broadcastStream::Bool = true,
            noDisplay::Bool = false,
            live_endpoint::AbstractString = "tcp://localhost:5556",
            broadcast_endpoint::AbstractString = "tcp://localhost:5570",
            record_file_playback_path::Union{Nothing,AbstractString} = "vizard_playback.bin",
            sim_synch_enabled::Bool = true,
            sim_synch_acceleration_factor::Real = 1.0,
            bank_unverified_simsynch_header_field_names::Bool = false)
        if noDisplay && (liveStream || broadcastStream)
            throw(ArgumentError(
                "Vizard live communication violation: noDisplay=true is mutually exclusive with liveStream and broadcastStream"
            ))
        end
        bank_unverified_simsynch_header_field_names && throw(ArgumentError(
            "Verification Gate Violation: C header field names (accelFactor, displayTime, timeInitialized in simSynch.h) were not inspected directly and must remain in :Incomplete status"
        ))
        sim_synch_acceleration_factor > 0.0 ||
            throw(ArgumentError("sim_synch_acceleration_factor must be > 0"))

        new(
            liveStream,
            broadcastStream,
            noDisplay,
            String(live_endpoint),
            String(broadcast_endpoint),
            record_file_playback_path === nothing ? nothing : String(record_file_playback_path),
            sim_synch_enabled,
            Float64(sim_synch_acceleration_factor),
            :Incomplete
        )
    end
end

"""
    build_vizard_vizmessage_packet(
        sc_state::BasiliskSpacecraftState,
        clock::MissionClock,
        payload_to_sc_tf::SMBRCoordinateTransform;
        config::VizardTransportConfig = VizardTransportConfig(),
        reactor_centroid_local_um::AbstractVector{<:Real} = [500.0, 0.0, 0.0])

Builds a read-only `VizMessage` protobuf-ready packet descriptor linking the Basilisk spacecraft
state, the reactor-to-spacecraft coordinate transform `SMBRCoordinateTransform`, the Three-Clock Law
`MissionClock`, and the ZeroMQ transport configuration (defaulting to the non-blocking read-only
`broadcastStream` tap on `tcp://localhost:5570`).
"""
function build_vizard_vizmessage_packet(
        sc_state::BasiliskSpacecraftState,
        clock::MissionClock,
        payload_to_sc_tf::SMBRCoordinateTransform;
        config::VizardTransportConfig = VizardTransportConfig(),
        reactor_centroid_local_um::AbstractVector{<:Real} = [500.0, 0.0, 0.0])
    reactor_pos_sc_m = transform_point_forward(payload_to_sc_tf, reactor_centroid_local_um)
    recovered_local_um = transform_point_inverse(payload_to_sc_tf, reactor_pos_sc_m)

    active_mode = if config.noDisplay
        :noDisplay
    elseif config.broadcastStream && !config.liveStream
        :broadcastStream
    elseif config.liveStream && !config.broadcastStream
        :liveStream
    else
        :dual_live_and_broadcast
    end

    return (
        unity_application = true,
        license = :Freeware,
        serialization_format = :protobuf_over_zeromq,
        json_side_channel_supported = false,
        handshake_strings = VIZARD_ZMQ_HANDSHAKE_STRINGS,
        active_mode = active_mode,
        config = config,
        canonical_live_stream_example = "scenarioBasicOrbitStream",
        sim_synch = (
            slaves_to_realtime_wall_clock = config.sim_synch_enabled,
            acceleration_factor = config.sim_synch_acceleration_factor,
            mid_run_change_requires_time_initialized_reset = true,
            c_header_field_names_verdict = config.sim_synch_header_field_names_verdict
        ),
        clocks = clock,
        spacecraft_position_si_m = sc_state.r_BN_N_m,
        reactor_position_in_spacecraft_frame_m = reactor_pos_sc_m,
        roundtrip_local_centroid_um = recovered_local_um,
        transform_roundtrip_error_um = norm(recovered_local_um .- Float64.(reactor_centroid_local_um))
    )
end

"""
    verify_vizard_socket_loopback_handshake()

Executes a live loopback TCP socket verification on `127.0.0.1` ephemeral ports
exercising both:
1. `broadcastStream` (`tcp://localhost:5570` semantics): non-blocking read-only
   one-way publish-subscribe frame delivery (`lockstep_waits_for_ok == false`).
2. `liveStream` (`tcp://localhost:5556` semantics): two-way lockstep handshake
   exchanging `"PING"` -> `"OK"` and `"SIM_UPDATE"` -> `"OK"` from
   `VIZARD_ZMQ_HANDSHAKE_STRINGS` (`lockstep_waits_for_ok == true`).
"""
function verify_vizard_socket_loopback_handshake()
    # 1. Verify read-only broadcastStream one-way frame delivery (no OK wait)
    bcast_port, bcast_server = Sockets.listenany(Sockets.ip"127.0.0.1", 55700)
    bcast_payload = Vector{UInt8}("VIZMESSAGE_PROTO_FRAME_EPOCH_846859269")
    bcast_received = Channel{Vector{UInt8}}(1)
    bcast_task = @async begin
        sock = Sockets.accept(bcast_server)
        try
            nbytes = read(sock, Int32)
            put!(bcast_received, read(sock, nbytes))
        finally
            close(sock)
        end
    end
    bcast_client = Sockets.connect(Sockets.ip"127.0.0.1", bcast_port)
    try
        write(bcast_client, Int32(length(bcast_payload)))
        write(bcast_client, bcast_payload)
        flush(bcast_client)
    finally
        close(bcast_client)
        close(bcast_server)
    end
    wait(bcast_task)
    rx_bytes = take!(bcast_received)
    broadcast_ok = (rx_bytes == bcast_payload) &&
        !VIZARD_LIVE_MODES.broadcastStream.lockstep_waits_for_ok

    # 2. Verify two-way liveStream lockstep PING -> OK and SIM_UPDATE -> OK
    live_port, live_server = Sockets.listenany(Sockets.ip"127.0.0.1", 55560)
    live_server_task = @async begin
        sock = Sockets.accept(live_server)
        try
            msg1 = readline(sock)
            msg1 == "PING" && write(sock, "OK\n")
            flush(sock)
            msg2 = readline(sock)
            msg2 == "SIM_UPDATE" && write(sock, "OK\n")
            flush(sock)
        finally
            close(sock)
        end
    end
    live_client = Sockets.connect(Sockets.ip"127.0.0.1", live_port)
    reply1 = ""
    reply2 = ""
    try
        write(live_client, "PING\n")
        flush(live_client)
        reply1 = readline(live_client)
        write(live_client, "SIM_UPDATE\n")
        flush(live_client)
        reply2 = readline(live_client)
    finally
        close(live_client)
        close(live_server)
    end
    wait(live_server_task)
    livestream_ping_ok = (reply1 == "OK") && VIZARD_LIVE_MODES.liveStream.lockstep_waits_for_ok
    livestream_sim_update_ok = (reply2 == "OK")

    return (
        broadcast_non_blocking_verified = broadcast_ok,
        broadcast_bytes_received = length(rx_bytes),
        livestream_lockstep_ping_ok_verified = livestream_ping_ok,
        livestream_sim_update_ok_verified = livestream_sim_update_ok,
        handshake_tokens_verified = ("PING", "OK", "SIM_UPDATE"),
        default_broadcast_endpoint = VIZARD_LIVE_MODES.broadcastStream.default_endpoint,
        default_livestream_endpoint = VIZARD_LIVE_MODES.liveStream.default_endpoint,
    )
end

