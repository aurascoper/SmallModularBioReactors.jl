@testset "Integration, Visualization, Verification Gates, and Milestones S0..S6" begin
    # 1. Conservation & NeuroARC 3-Valued Gates
    cons = audit_reactor_conservation_invariants()
    @test cons.all_conserved
    @test cons.role_contract.assigned_role == :evidence_verifier

    g_acc = evaluate_neuroarc_evidence_gate(claim_id = "witnessed_claim", has_witness_evidence = true, controls_satisfied = true)
    g_inc = evaluate_neuroarc_evidence_gate(claim_id = "unwitnessed_claim", has_witness_evidence = false)
    g_rej = evaluate_neuroarc_evidence_gate(claim_id = "refuted_claim", has_refutation_evidence = true)
    @test g_acc.verdict == :Accepted && g_acc.eligible_for_memory_banking
    @test g_inc.verdict == :Incomplete && !g_inc.eligible_for_memory_banking
    @test g_rej.verdict == :Rejected && !g_rej.eligible_for_memory_banking
    @test g_inc.incomplete_never_collapses_to_rejected
    @test_throws ArgumentError evaluate_neuroarc_evidence_gate(
        claim_id = "illegal_collapse",
        has_witness_evidence = false,
        force_collapse_incomplete_to_rejected = true
    )

    memo_audit = audit_verified_memo_memory_banking_rules()
    @test memo_audit.all_rules_enforced
    @test_throws ArgumentError audit_verified_memo_memory_banking_rules(bank_unconfirmed_claim = true)

    # 2. Uncertainty & Equifinality
    unc = decompose_epistemic_aleatoric_uncertainty([1.0 1.2 0.8; 2.0 2.1 1.9], 0.01)
    @test all(unc.total_predictive_variance .> unc.aleatoric_variance)
    eqf = audit_equifinality_underdetermination([0.6, 0.4], [0.3, 0.7], x -> [sum(x)])
    @test eqf.is_equifinal && eqf.unique_backward_inversion_forbidden

    # 3. Visualization: M7Adapter & VizardAdapter (3 Live Modes & Handshake Strings)
    vcfg_bcast = VizardTransportConfig(liveStream = false, broadcastStream = true, noDisplay = false)
    @test vcfg_bcast.broadcast_endpoint == "tcp://localhost:5570"
    @test vcfg_bcast.sim_synch_header_field_names_verdict == :Incomplete
    @test_throws ArgumentError VizardTransportConfig(liveStream = true, noDisplay = true)
    @test_throws ArgumentError VizardTransportConfig(broadcastStream = true, noDisplay = true)
    @test_throws ArgumentError VizardTransportConfig(bank_unverified_simsynch_header_field_names = true)

    pkt = build_vizard_vizmessage_packet(BasiliskSpacecraftState(), MissionClock(), SMBRCoordinateTransform(); config = vcfg_bcast)
    @test pkt.unity_application
    @test pkt.serialization_format == :protobuf_over_zeromq
    @test !pkt.json_side_channel_supported
    @test pkt.handshake_strings == ("PING", "SIM_UPDATE", "OK", "REQUEST_INPUT", "VIZARD_INPUT", "REQUEST_IMAGE_")
    @test pkt.active_mode == :broadcastStream
    @test pkt.transform_roundtrip_error_um < 1.0e-9

    # 4. Demonstrations 1..3 and Milestones S0..S6
    suite = evaluate_smbr_milestones_s0_to_s6()
    @test suite.all_milestones_passed
    @test suite.demonstrations.d1.passed
    @test suite.demonstrations.d2.passed
    @test suite.demonstrations.d3.passed

    # 5. Strictly zero em-dashes and zero en-dashes across the entire package
    pkg_root = normpath(joinpath(@__DIR__, "..", ".."))
    for (root, _, files) in walkdir(pkg_root)
        for f in files
            if endswith(f, ".jl") || endswith(f, ".md") || endswith(f, ".toml")
                txt = read(joinpath(root, f), String)
                @test !occursin("\u2014", txt)
                @test !occursin("\u2013", txt)
            end
        end
    end
end
