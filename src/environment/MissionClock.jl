# ==============================================================================
# src/environment/MissionClock.jl
# Three-Clock Law: (1) Mission Elapsed Time (SPICE/Basilisk), (2) Physical Flow
# Time (OpenFOAM/DBS), (3) Geochemical/Biological Kinetic Time (PHREEQC/Biofilm),
# plus Discrete Sobolev Relaxation Step
# ==============================================================================

struct MissionClock
    mission_elapsed_time_s::Float64
    physical_flow_time_s::Float64
    kinetic_biogeochem_time_s::Float64
    sobolev_relaxation_step::Int
    sync_event_id::Union{Nothing,String}
    clocks_explicitly_coupled::Bool

    function MissionClock(;
            mission_elapsed_time_s::Real = 0.0,
            physical_flow_time_s::Real = 0.0,
            kinetic_biogeochem_time_s::Real = 0.0,
            sobolev_relaxation_step::Int = 0,
            sync_event_id::Union{Nothing,AbstractString} = nothing,
            clocks_explicitly_coupled::Bool = false)
        (isfinite(mission_elapsed_time_s) && mission_elapsed_time_s >= 0.0) ||
            throw(ArgumentError("mission_elapsed_time_s must be finite and >= 0"))
        (isfinite(physical_flow_time_s) && physical_flow_time_s >= 0.0) ||
            throw(ArgumentError("physical_flow_time_s must be finite and >= 0"))
        (isfinite(kinetic_biogeochem_time_s) && kinetic_biogeochem_time_s >= 0.0) ||
            throw(ArgumentError("kinetic_biogeochem_time_s must be finite and >= 0"))
        sobolev_relaxation_step >= 0 ||
            throw(ArgumentError("sobolev_relaxation_step must be >= 0"))
        if clocks_explicitly_coupled && (sync_event_id === nothing || isempty(strip(sync_event_id)))
            throw(ArgumentError("Three-Clock Law violation: clocks_explicitly_coupled=true requires a non-empty sync_event_id"))
        end
        new(
            Float64(mission_elapsed_time_s),
            Float64(physical_flow_time_s),
            Float64(kinetic_biogeochem_time_s),
            sobolev_relaxation_step,
            sync_event_id === nothing ? nothing : String(sync_event_id),
            clocks_explicitly_coupled
        )
    end
end

"""
    advance_mission_clock(clock::MissionClock;
        d_mission_s::Real = 0.0,
        d_flow_s::Real = 0.0,
        d_kinetic_s::Real = 0.0,
        d_sobolev_steps::Int = 0,
        sync_event_id::Union{Nothing,AbstractString} = nothing,
        treat_sobolev_as_physical_time::Bool = false,
        assume_clocks_equal_without_event::Bool = false)

Advances the three continuous clocks and discrete Sobolev step independently while enforcing
the Three-Clock Law. Throws `ArgumentError` if Sobolev relaxation steps are treated as physical
seconds or if clocks are assumed equal without an explicit synchronization event.
"""
function advance_mission_clock(
        clock::MissionClock;
        d_mission_s::Real = 0.0,
        d_flow_s::Real = 0.0,
        d_kinetic_s::Real = 0.0,
        d_sobolev_steps::Int = 0,
        sync_event_id::Union{Nothing,AbstractString} = nothing,
        treat_sobolev_as_physical_time::Bool = false,
        assume_clocks_equal_without_event::Bool = false)
    treat_sobolev_as_physical_time &&
        throw(ArgumentError("Three-Clock Law violation: Sobolev relaxation steps cannot be treated as physical seconds"))
    if assume_clocks_equal_without_event && (sync_event_id === nothing || isempty(strip(sync_event_id)))
        throw(ArgumentError("Three-Clock Law violation: clocks cannot be assumed equal without an explicit sync_event_id"))
    end
    (d_mission_s >= 0.0 && d_flow_s >= 0.0 && d_kinetic_s >= 0.0 && d_sobolev_steps >= 0) ||
        throw(ArgumentError("Clock increments must be non-negative"))

    coupled = (sync_event_id !== nothing && !isempty(strip(sync_event_id)))
    return MissionClock(
        mission_elapsed_time_s = clock.mission_elapsed_time_s + Float64(d_mission_s),
        physical_flow_time_s = clock.physical_flow_time_s + Float64(d_flow_s),
        kinetic_biogeochem_time_s = clock.kinetic_biogeochem_time_s + Float64(d_kinetic_s),
        sobolev_relaxation_step = clock.sobolev_relaxation_step + d_sobolev_steps,
        sync_event_id = coupled ? String(sync_event_id) : clock.sync_event_id,
        clocks_explicitly_coupled = coupled
    )
end
