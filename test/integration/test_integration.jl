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

    # 5. Health-Physics ALARA (Time, Distance 1/r^2, Shielding + Melanin) & Epistemic Hippocratism
    alara_nom = evaluate_health_physics_alara_and_hippocratism()
    @test alara_nom.alara_satisfied
    @test alara_nom.primum_non_nocere_enforced
    @test alara_nom.epistemic_hippocratism_verdict == :Accepted
    @test isapprox(alara_nom.distance_control.inverse_square_attenuation_factor, 0.01; atol = 1.0e-15)
    @test alara_nom.shielding_control.alara_dose_reduction_factor > 100.0
    @test !alara_nom.time_control.payload_safing_triggered

    alara_cold = evaluate_health_physics_alara_and_hippocratism(
        SpectralRadiationFieldSpec(),
        HealthPhysicsALARAContract();
        temperature_K = 270.0
    )
    @test alara_cold.below_280K_underconstrained
    @test alara_cold.epistemic_hippocratism_verdict == :Incomplete

    alara_unverified_spice = evaluate_health_physics_alara_and_hippocratism(
        SpectralRadiationFieldSpec(),
        HealthPhysicsALARAContract();
        spice_units_verified_by_spkpos = false
    )
    @test alara_unverified_spice.epistemic_hippocratism_verdict == :Incomplete

    @test_throws ArgumentError HealthPhysicsALARAContract(standoff_distance_m = 0.5, reference_distance_m = 1.0)
    @test_throws ArgumentError HealthPhysicsALARAContract(epistemic_hippocratism_enforced = false)
    @test_throws ArgumentError evaluate_health_physics_alara_and_hippocratism(claim_unverified_low_dose_hormesis_as_safe = true)
    @test_throws ArgumentError evaluate_health_physics_alara_and_hippocratism(recursion_depth = 1)

    # 6. M7.1..M7.8 Header vs Float64/Int64 Schema, Selfish-Gene/Comfort Reticulated ASR, and Tarjan Bridge DAG
    m7_sch = validate_m7_header_vs_numeric_schema()
    @test m7_sch.schema_valid
    @test m7_sch.stage_headers == VALID_M7_STAGE_HEADERS
    @test m7_sch.ieee754_binary_approximation_residual > 0.0
    @test_throws ArgumentError validate_m7_header_vs_numeric_schema(conflate_header_as_numeric_literal = true)
    @test_throws ArgumentError validate_m7_header_vs_numeric_schema(stage_headers = ["M7.1", "M7.2"])

    asr_ret = evaluate_selfish_gene_and_comfort_reticulated_asr(hgt_reticulation_edges = 2)
    @test asr_ret.reticulated_asr_verified
    @test asr_ret.beta_1_reticulation_cycles == 2
    @test asr_ret.is_reticulated_network
    @test asr_ret.inter_level_tree_disagreement_is_signal
    @test_throws ArgumentError evaluate_selfish_gene_and_comfort_reticulated_asr(
        hgt_reticulation_edges = 2,
        force_strict_tree_on_reticulated_network = true
    )

    t_dag = evaluate_tarjan_bridge_julia_matlab_pipeline()
    @test t_dag.is_acyclic_tarjan_bridge_chain
    @test t_dag.n_tarjan_bridge_edges == 4
    @test t_dag.beta_1_inter_solver == 0
    @test !t_dag.matlab_mutates_upstream_solver
    @test_throws ArgumentError evaluate_tarjan_bridge_julia_matlab_pipeline(introduce_feedback_cycle_to_openfoam = true)

    # 7. Strictly zero em-dashes and zero en-dashes across the entire package
    pkg_root = normpath(joinpath(@__DIR__, "..", ".."))
    for (root, _, files) in walkdir(pkg_root)
        for f in files
            if endswith(f, ".jl") || endswith(f, ".md") || endswith(f, ".toml") || endswith(f, ".m") || endswith(f, ".json")
                txt = read(joinpath(root, f), String)
                @test !occursin("\u2014", txt)
                @test !occursin("\u2013", txt)
            end
        end
    end
end
