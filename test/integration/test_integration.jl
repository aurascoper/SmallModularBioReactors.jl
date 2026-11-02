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
    @test !alara_nom.distance_control.at_exact_standoff_boundary
    @test alara_nom.shielding_control.alara_dose_reduction_factor > 100.0
    @test !alara_nom.time_control.payload_safing_triggered

    # Edge case: exact geometric standoff boundary r == r_0 and exact NutMEG lower boundary T == 280.0 K
    alara_exact_r0_280K = evaluate_health_physics_alara_and_hippocratism(
        SpectralRadiationFieldSpec(),
        HealthPhysicsALARAContract(reference_distance_m = 1.0, standoff_distance_m = 1.0);
        temperature_K = 280.0
    )
    @test alara_exact_r0_280K.distance_control.at_exact_standoff_boundary
    @test isapprox(alara_exact_r0_280K.distance_control.inverse_square_attenuation_factor, 1.0; atol = 1.0e-15)
    @test alara_exact_r0_280K.at_exact_280K_nutmeg_boundary
    @test !alara_exact_r0_280K.below_280K_underconstrained
    @test alara_exact_r0_280K.epistemic_hippocratism_verdict == :Accepted

    # Edge case: cumulative dose inside [0.80 * D_max, D_max] payload safing band + T == 280.0 K
    peak_rate = maximum(alara_nom.dose_rate_Gy_s)
    t_safing_band = (0.85 * 0.50) / peak_rate
    alara_safing_band = evaluate_health_physics_alara_and_hippocratism(
        SpectralRadiationFieldSpec(),
        HealthPhysicsALARAContract(time_exposure_s = t_safing_band, max_allowable_cumulative_dose_Gy = 0.50);
        temperature_K = 280.0
    )
    @test alara_safing_band.time_control.payload_safing_triggered
    @test !alara_safing_band.alara_satisfied
    @test alara_safing_band.epistemic_hippocratism_verdict == :Incomplete

    # Edge case: cumulative dose exceeding D_max -> :Rejected
    t_exceed_band = (1.15 * 0.50) / peak_rate
    alara_exceeded = evaluate_health_physics_alara_and_hippocratism(
        SpectralRadiationFieldSpec(),
        HealthPhysicsALARAContract(time_exposure_s = t_exceed_band, max_allowable_cumulative_dose_Gy = 0.50);
        temperature_K = 298.15
    )
    @test !alara_exceeded.alara_satisfied
    @test alara_exceeded.epistemic_hippocratism_verdict == :Rejected

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
    @test_throws ArgumentError evaluate_health_physics_alara_and_hippocratism(depths_m = Float64[])

    # 5b. PayloadController, MelaninResponse ALARA Shielding, and EvidenceGates Epistemic Hippocratism
    pc_nom = evaluate_alara_payload_controller_safing()
    @test pc_nom.controller_mode == :NOMINAL_BIOREACTOR_SCIENCE
    @test !pc_nom.shutter_closed
    @test pc_nom.alara_satisfied
    @test pc_nom.epistemic_verdict == :Accepted

    pc_safe = evaluate_alara_payload_controller_safing(
        MissionClock(mission_elapsed_time_s = t_safing_band),
        HealthPhysicsALARAContract(time_exposure_s = t_safing_band, max_allowable_cumulative_dose_Gy = 0.50)
    )
    @test pc_safe.controller_mode == :SAFE_HOLD_SHUTTER_CLOSED
    @test pc_safe.shutter_closed
    @test pc_safe.epistemic_verdict == :Incomplete

    pc_scram = evaluate_alara_payload_controller_safing(
        MissionClock(mission_elapsed_time_s = t_exceed_band),
        HealthPhysicsALARAContract(time_exposure_s = t_exceed_band, max_allowable_cumulative_dose_Gy = 0.50)
    )
    @test pc_scram.controller_mode == :EMERGENCY_SCRAM_ISOLATE
    @test pc_scram.shutter_closed
    @test pc_scram.epistemic_verdict == :Rejected

    @test_throws ArgumentError evaluate_alara_payload_controller_safing(algorithmic_exit_trigger_stoppable = false)
    @test_throws ArgumentError PayloadSafetyControllerSpec(automated_shutter_isolation_enabled = false)

    mel_alara = evaluate_melanin_neurospora_alara_shielding_factor()
    @test mel_alara.alara_biological_protection_verified
    @test mel_alara.radical_suppression_ratio > 1.5
    @test_throws ArgumentError evaluate_melanin_neurospora_alara_shielding_factor(relax_physical_shielding_via_hormesis = true)

    eg_hippo = evaluate_epistemic_hippocratism_quarantine_gate()
    @test eg_hippo.verdict == :Accepted
    @test eg_hippo.primum_non_nocere_enforced
    @test eg_hippo.eligible_for_memory_banking
    eg_hippo_safing = evaluate_epistemic_hippocratism_quarantine_gate(cumulative_dose_Gy = 0.42, max_allowable_dose_Gy = 0.50)
    @test eg_hippo_safing.verdict == :Incomplete
    @test eg_hippo_safing.safing_band_quarantined
    @test !eg_hippo_safing.eligible_for_memory_banking
    @test_throws ArgumentError evaluate_epistemic_hippocratism_quarantine_gate(promote_unverified_hormesis_or_radiosynthesis = true)
    @test_throws ArgumentError evaluate_epistemic_hippocratism_quarantine_gate(recursion_depth = 2)

    # 6. M7.1..M7.8 Header vs Float64/Int64 Schema, Selfish-Gene/Comfort Reticulated ASR, Tarjan Bridge DAG, and MATLAB Reactor Pools / Biofilms Seam
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

    openmc_l32 = evaluate_openmc_cross_section_substrate_swap_gate()
    @test openmc_l32.valid
    @test openmc_l32.library_primary == "ENDF/B-VIII.0"
    @test openmc_l32.library_counterfactual == "JEFF-3.3"
    @test openmc_l32.tally_variance_rel_diff <= 0.05
    @test !evaluate_openmc_cross_section_substrate_swap_gate(dose_rate_endfb80_mSv_hr = 0.010, dose_rate_jeff33_mSv_hr = 0.020).valid
    @test_throws ArgumentError evaluate_openmc_cross_section_substrate_swap_gate(icrp116_dosimetry_verified_separately = false)
    @test_throws ArgumentError evaluate_openmc_cross_section_substrate_swap_gate(workplace_regulatory_geometry_verified = false)
    @test_throws ArgumentError evaluate_openmc_cross_section_substrate_swap_gate(recursion_depth = 1)

    polymesh_rt = parse_openfoam_polymesh_topology_and_roundtrip()
    @test polymesh_rt.roundtrip_beta1_invariant
    @test polymesh_rt.n_cells === Int64(4)
    @test polymesh_rt.n_internal_faces === Int64(3)
    @test polymesh_rt.n_cyclic_collar_patches === Int64(2)
    @test polymesh_rt.beta_0_openfoam === Int64(1)
    @test polymesh_rt.beta_1_openfoam === Int64(2)
    @test polymesh_rt.field_interpolation_tear_free
    @test_throws ArgumentError parse_openfoam_polymesh_topology_and_roundtrip(expected_beta_1 = 99)
    @test_throws ArgumentError parse_openfoam_polymesh_topology_and_roundtrip(matlab_sampled_field = [1.0, 2.0])

    mp_seam = evaluate_matlab_neuralcompose_collar_and_multipool_biofilm_seam()
    @test mp_seam.collar_and_multipool_valid
    @test mp_seam.explicit_cross_process_collar_enforced
    @test mp_seam.mat_v7_binary_supported && mp_seam.mat_v73_hdf5_supported
    @test mp_seam.mat_v7_magic_prefix == "MATLAB 5.0 MAT-file"
    @test mp_seam.mat_v73_hdf5_magic_hex == "894844460d0a1a0a"
    @test mp_seam.openfoam_basilisk_matlab_roundtrip_beta1_invariant
    @test mp_seam.openfoam_polymesh_roundtrip.beta_1_openfoam === Int64(2)
    @test mp_seam.headless_cli_diff_checker_verified
    @test mp_seam.swarm_biconnected_verification_valid
    @test mp_seam.beta_1_swarm_verification_cycles == 2
    @test mp_seam.counterfactual_feed_swap_trigger_invariant
    @test mp_seam.unified_substrate_swap_L31_valid
    @test mp_seam.stage_rejects_logged_separately_L31
    @test mp_seam.openmc_cross_section_swap_L32.valid
    @test mp_seam.openmc_tally_variance_rel_diff <= 0.05
    @test mp_seam.multica_synthetic_warrant_models == ("claude-sonnet-5", "claude-sonnet-5-5", "claude-haiku-4-5")
    @test mp_seam.multipool_biofilm_seam.biomass_conserved
    @test mp_seam.multipool_biofilm_seam.n_pools == 4
    @test mp_seam.n_modular_chambers === Int64(8)
    @test mp_seam.multipool_8chamber_seam.n_pools == 8
    @test mp_seam.multipool_8chamber_seam.biomass_conserved
    @test mp_seam.beta_0_swarm_components === Int64(1)
    @test isempty(mp_seam.swarm_cut_vertices)
    @test isempty(mp_seam.swarm_bridge_edges)
    @test !isempty(mp_seam.stage1_mesh_rejects_L31)
    @test !isempty(mp_seam.stage2_swarm_rejects_L31)
    @test !isempty(mp_seam.stage3_openmc_rejects_L31)
    @test mp_seam.neurospora_8pool_state.all_4_blocks_valid
    @test_throws ArgumentError evaluate_matlab_neuralcompose_collar_and_multipool_biofilm_seam(assume_shared_in_memory_process_without_collar = true)
    @test_throws ArgumentError evaluate_matlab_neuralcompose_collar_and_multipool_biofilm_seam(openfoam_matlab_roundtrip_beta1_preserved = false)
    @test_throws ArgumentError evaluate_matlab_neuralcompose_collar_and_multipool_biofilm_seam(swarm_has_single_point_bridge_failure = true)
    @test_throws ArgumentError evaluate_matlab_neuralcompose_collar_and_multipool_biofilm_seam(stage_rejects_logged_separately_L31 = false)
    @test_throws ArgumentError evaluate_matlab_neuralcompose_collar_and_multipool_biofilm_seam(mat_v73_header_bytes = Vector{UInt8}("CORRUPTED_NOT_HDF5"))
    @test_throws ArgumentError evaluate_matlab_neuralcompose_collar_and_multipool_biofilm_seam(
        swarm_nodes = [:a, :b, :c, :d, :e],
        swarm_edges = [:a => :b, :b => :c, :c => :a, :c => :d, :d => :e, :e => :c]
    )
    @test_throws ArgumentError coordinate_matlab_multipool_biofilm_session_seam(overwrite_peer_session_state = true)

    nc8_preset = evaluate_neurospora_crassa_syncytial_circadian_preset(
        NeurosporaCrassaCircadianHyphalPreset(n_modular_chambers = 8)
    )
    @test nc8_preset.syncytial_signal_propagation_intact
    @test length(nc8_preset.chamber_hyphal_extension_um_hr) == 8

    nc4_state = evaluate_neurospora_crassa_4block_matlab_state()
    @test nc4_state.all_4_blocks_valid
    @test nc4_state.n_pool === Int64(8)
    @test length(nc4_state.M_frq) == 8 && length(nc4_state.F_c) == 8 && length(nc4_state.F_n) == 8 && length(nc4_state.phi_t) == 8
    @test length(nc4_state.B_i) == 8 && length(nc4_state.rho_tip_i) == 8 && length(nc4_state.phi_EPS_i) == 8
    @test nc4_state.beta_0 === Int64(1) && nc4_state.beta_1 === Int64(3)
    @test length(nc4_state.S_i) == 8 && length(nc4_state.C_O2_i) == 8 && length(nc4_state.D_mel_i) == 8 && length(nc4_state.D_dot_eff_i) == 8
    @test length(nc4_state.conidial_cytometry.FSC_A) == 8 && length(nc4_state.conidial_cytometry.SSC_A) == 8
    @test length(nc4_state.conidial_cytometry.FL1_A) == 8 && length(nc4_state.conidial_cytometry.melanin_pg_cell) == 8
    @test nc4_state.openmc_tally_variance_rel_diff <= 0.05
    @test length(nc4_state.u_valve) == 8 && length(nc4_state.z_shield) == 8 && length(nc4_state.Phi_light) == 8
    @test nc4_state.nutmeg_verdict == :Accepted
    @test nc4_state.hormesis_quarantine_verdict == :Incomplete
    @test evaluate_neurospora_crassa_4block_matlab_state(temperature_K = 275.0).nutmeg_verdict == :Incomplete
    @test_throws ArgumentError evaluate_neurospora_crassa_4block_matlab_state(promote_radiation_hormesis = true)

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
