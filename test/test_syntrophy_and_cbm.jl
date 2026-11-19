@testset "Syntrophic Bioenergetics (SYN0..SYN6) & M7.9 Clay-Biofilm-Metal (CBM-1)" begin
    # ------------------------------------------------------------------
    # SYN0: Stoichiometric & electron-equivalent ledger
    # ------------------------------------------------------------------
    s0_prop = evaluate_syntrophic_stoichiometric_ledger(
        substrate = :propionate,
        donor_consumed_mol = 0.0203
    )
    @test s0_prop.electron_balance_closed
    @test s0_prop.electron_balance_rel_residual < 1.0e-12
    @test isapprox(s0_prop.ch4_produced_mol, 0.015225; atol = 1.0e-12)
    @test !s0_prop.butanol_shift_active

    # Alternative acceptors (SO4^2-, NO3-) + formate + >390 ppm H2 butanol shift
    s0_alt = evaluate_syntrophic_stoichiometric_ledger(
        substrate = :butyrate,
        donor_consumed_mol = 0.020,
        h2_accumulated_mol = 0.002,
        formate_accumulated_mol = 0.001,
        alt_acceptor_so4_reduced_mol = 0.001,
        alt_acceptor_no3_reduced_mol = 0.002,
        p_h2_ppm = 450.0,
        bicarbonate_limited = true,
        butanol_produced_mol = 0.003
    )
    @test s0_alt.electron_balance_closed
    @test s0_alt.electron_balance_rel_residual < 1.0e-12
    @test s0_alt.butanol_shift_triggered_above_390ppm
    @test s0_alt.butanol_shift_active
    @test s0_alt.e_alternative_acceptors_mol > 0.0

    # Unbalanced electron ledger & bogus non-tautological carbon/acetate ledger refused/detected
    s0_bad = evaluate_syntrophic_stoichiometric_ledger(
        substrate = :propionate,
        donor_consumed_mol = 0.020,
        ch4_produced_mol = 0.005
    )
    @test !s0_bad.electron_balance_closed
    s0_bogus_ac = evaluate_syntrophic_stoichiometric_ledger(
        substrate = :propionate,
        donor_consumed_mol = 0.020,
        acetate_produced_mol = 999.0
    )
    @test !s0_bogus_ac.carbon_balance_closed
    @test !s0_bogus_ac.acetate_balance_closed
    @test !s0_bogus_ac.electron_balance_closed
    @test s0_bogus_ac.carbon_balance_rel_residual > 0.9
    @test_throws ArgumentError evaluate_syntrophic_stoichiometric_ledger(
        substrate = :propionate,
        donor_consumed_mol = 0.020,
        butanol_produced_mol = 0.015
    )
    @test_throws ArgumentError evaluate_syntrophic_stoichiometric_ledger(force_unbalanced_electron_ledger = true)

    # ------------------------------------------------------------------
    # SYN1: Partner qPCR enforcement + Monod vs Moser PHF-HUR (70 C cutoff)
    # ------------------------------------------------------------------
    @test_throws ArgumentError step_syntrophic_two_population_kinetics(use_total_biomass_only = true)
    s1_37 = step_syntrophic_two_population_kinetics(temperature_C = 37.0, kinetic_model = :moser)
    s1_65 = step_syntrophic_two_population_kinetics(temperature_C = 65.0, kinetic_model = :moser)
    s1_70 = step_syntrophic_two_population_kinetics(temperature_C = 70.0, kinetic_model = :moser)

    @test s1_37.partner_qpcr_enforced
    @test s1_37.electron_balance_conserved
    @test s1_37.moser_preferred_over_monod
    @test s1_37.q_e_methanogen_mol_kg_s > s1_37.q_e_donor_mol_kg_s
    @test isapprox(s1_65.q_max_H2_mmol_gVSS_d, 773.0; rtol = 1.0e-6)
    @test s1_70.thermal_cutoff_at_70C
    @test s1_70.q_max_H2_mol_kg_s == 0.0
    @test s1_70.moser_specific_rate_mol_kg_s == 0.0
    @test s1_70.CH4_produced_mol_m3 == 0.0

    # ------------------------------------------------------------------
    # SYN2: Gas/liquid H2 transfer & k_L*a vs PHF-HUR K_S bias
    # ------------------------------------------------------------------
    s2_gas = evaluate_syntrophic_gas_liquid_transfer(headspace_free_phf_hur = false)
    s2_phf = evaluate_syntrophic_gas_liquid_transfer(headspace_free_phf_hur = true)
    @test s2_gas.gas_liquid_transfer_limited
    @test s2_gas.ks_bias_ratio > 1.0
    @test s2_gas.K_S_apparent_gas_supplied_mol_m3 > s2_gas.K_S_true_mol_m3
    @test s2_gas.liquid_oversaturated
    @test s2_gas.K_S_apparent_liquid_oversat_mol_m3 < s2_gas.K_S_true_mol_m3
    @test s2_phf.phf_hur_eliminates_kla_bias
    @test isapprox(s2_phf.ks_bias_ratio, 1.0; atol = 1.0e-12)

    # ------------------------------------------------------------------
    # SYN3: Thermodynamics, NutMEG closure, & regulatory silence separation
    # ------------------------------------------------------------------
    s3_ok = evaluate_syntrophic_thermodynamics(temperature_K = 310.15, p_H2_atm = 1.0e-5)
    @test s3_ok.donor_oxidation_exergonic
    @test s3_ok.methanogenesis_exergonic
    @test s3_ok.syntrophy_thermodynamically_feasible
    @test !s3_ok.thermodynamic_feasibility_sufficient_for_monoculture_h2
    @test s3_ok.verdict == :Accepted

    # High H2 shuts off propionate oxidation thermodynamically -> :Incomplete
    s3_high_h2 = evaluate_syntrophic_thermodynamics(temperature_K = 310.15, p_H2_atm = 0.1)
    @test !s3_high_h2.donor_oxidation_exergonic
    @test s3_high_h2.verdict == :Incomplete

    # Below 280 K NutMEG gate holds at :Incomplete
    s3_cold = evaluate_syntrophic_thermodynamics(temperature_K = 275.0, p_H2_atm = 1.0e-5)
    @test s3_cold.nutmeg_closure.out_of_calibrated_maintenance_range_below_280K
    @test s3_cold.verdict == :Incomplete

    @test_throws ArgumentError evaluate_syntrophic_thermodynamics(
        equate_thermodynamic_feasibility_with_monoculture_flux = true
    )

    # ------------------------------------------------------------------
    # SYN4: Fisher identifiability rank check & evolutionary quarantine
    # ------------------------------------------------------------------
    s4_full = evaluate_syntrophic_identifiability_and_calibration(
        has_qpcr_partner_biomass = true,
        has_dissolved_h2_probe = true,
        headspace_free_phf_hur = true
    )
    @test s4_full.identifiable
    @test s4_full.fisher_rank == s4_full.n_parameters
    @test s4_full.identifiability_verdict == :Accepted
    @test s4_full.evolutionary_track_separated
    @test s4_full.evolutionary_quarantine.hydrogen_hypothesis_ancestral_engulfment.verdict == :Incomplete
    @test s4_full.evolutionary_quarantine.asgard_heimdallarchaeia_obligate_h2_syntrophy.verdict == :Incomplete

    s4_no_qpcr = evaluate_syntrophic_identifiability_and_calibration(
        has_qpcr_partner_biomass = false,
        has_dissolved_h2_probe = true,
        headspace_free_phf_hur = true
    )
    @test !s4_no_qpcr.identifiable
    @test s4_no_qpcr.fisher_rank < s4_no_qpcr.n_parameters
    @test s4_no_qpcr.identifiability_verdict == :Incomplete

    s4_no_h2 = evaluate_syntrophic_identifiability_and_calibration(
        has_qpcr_partner_biomass = true,
        has_dissolved_h2_probe = false,
        headspace_free_phf_hur = false
    )
    @test !s4_no_h2.identifiable
    @test s4_no_h2.identifiability_verdict == :Incomplete
    @test_throws ArgumentError evaluate_syntrophic_identifiability_and_calibration(
        promote_evolutionary_narrative_to_accepted = true
    )

    # ------------------------------------------------------------------
    # SYN5: Spatial aggregate diffusion & CPM lattice refusal
    # ------------------------------------------------------------------
    s5_disp = evaluate_syntrophic_spatial_diffusion(interspecies_distance_um = 5.30)
    s5_coag = evaluate_syntrophic_spatial_diffusion(interspecies_distance_um = 0.29)
    @test s5_coag.spatial_coaggregation_increases_flux
    @test isapprox(s5_coag.fickian_flux_amplification_vs_5p30_um, 5.30 / 0.29; rtol = 1.0e-12)
    @test_throws ArgumentError evaluate_syntrophic_spatial_diffusion(
        interspecies_distance_um = 2.0,
        coordinate_unit = :cpm_lattice_sites,
        lattice_pitch_um = nothing
    )
    s5_cal_cpm = evaluate_syntrophic_spatial_diffusion(
        interspecies_distance_um = 2.0,
        coordinate_unit = :cpm_lattice_sites,
        lattice_pitch_um = 0.5
    )
    @test isapprox(s5_cal_cpm.interspecies_distance_um, 1.0; atol = 1.0e-12)

    # ------------------------------------------------------------------
    # M7.9-CBM: 5-field state, 5-pool S ledger, 4-term density, hydration, gates
    # ------------------------------------------------------------------
    cbm = ClayBiofilmMetalState()
    @test cbm.cell_count == 5
    @test all(0.0 < x < 1.0 for x in cbm.theta_mobile)
    @test all(0.0 < x <= 1.0 for x in cbm.aw)
    @test_throws ArgumentError ClayBiofilmMetalState(theta_mobile = [1.2, 0.3, 0.3, 0.3, 0.3])
    @test_throws ArgumentError ClayBiofilmMetalState(aw = [1.05, 0.98, 0.97, 0.95, 0.93])

    sled = evaluate_cbm_five_pool_sulfur_ledger()
    @test length(sled.pool_names) == 5
    @test sled.sulfur_conserved
    @test sled.mass_balance_rel_residual < 1.0e-12
    @test sled.iron_mineral_retention_distinct_from_system_removal
    @test_throws ArgumentError evaluate_cbm_five_pool_sulfur_ledger(
        equate_fe_mineral_retention_with_system_removal = true
    )

    dens = evaluate_clay_density_four_terms()
    @test dens.four_terms_kept_separate
    @test dens.skb_tr1609_activity_cutoff_wet_kg_m3.range_wet_kg_m3 == (1740.0, 1880.0)
    @test dens.skb_design_target_wet_kg_m3.range_wet_kg_m3 == (1950.0, 2050.0)
    @test !dens.skb_design_target_wet_kg_m3.is_biological_threshold
    @test dens.unresolved_secondary_triad.verdict == :Incomplete
    @test_throws ArgumentError evaluate_clay_density_four_terms(collapse_density_terms_to_single_threshold = true)

    h0w = evaluate_clay_hydration_and_ectem(hydration_state = :state_0W)
    h1w = evaluate_clay_hydration_and_ectem(hydration_state = :state_1W, cation = :Na_plus)
    h2w = evaluate_clay_hydration_and_ectem(hydration_state = :state_2W, cation = :Ca2_plus)
    h3w = evaluate_clay_hydration_and_ectem(hydration_state = :state_3W, cation = :Na_plus)
    hectem = evaluate_clay_hydration_and_ectem(hydration_state = :state_2W, fe3_bioreduced_fraction = 0.26, imaging_mode = :ec_tem_hydrated)
    hvac = evaluate_clay_hydration_and_ectem(hydration_state = :state_2W, imaging_mode = :conventional_vacuum_tem)
    h0w_red = evaluate_clay_hydration_and_ectem(hydration_state = :state_0W, fe3_bioreduced_fraction = 0.01, imaging_mode = :ec_tem_hydrated)
    @test h0w_red.observed_d001_nm <= h0w.baseline_d001_A / 10.0
    @test h0w.baseline_d001_A < h1w.baseline_d001_A < h2w.baseline_d001_A < h3w.baseline_d001_A
    @test h3w.osmotic_swelling_possible
    @test isapprox(hectem.ectem_d001_nm, 1.1; atol = 1.0e-12)
    @test hectem.verdict == :Accepted
    @test hvac.vacuum_tem_quarantined
    @test hvac.verdict == :Incomplete
    @test_throws ArgumentError evaluate_clay_hydration_and_ectem(promote_vacuum_tem_d001_as_hydrated = true)

    # Gates CBM-G1, CBM-G2, CBM-G3
    g1_cerna = evaluate_cbm_g1_abiotic_compaction_gate(
        dry_density_kg_m3 = 1600.0,
        abiotic_sterile_rate_um_yr = 2.4,
        biotic_rate_um_yr = 1.8
    )
    @test g1_cerna.abiotic_dominated
    @test g1_cerna.biotic_mic_attribution_verdict == :Rejected
    g1_below_lod = evaluate_cbm_g1_abiotic_compaction_gate(
        abiotic_sterile_rate_um_yr = 0.0,
        biotic_rate_um_yr = 0.0,
        detection_limit_um_yr = 0.25
    )
    @test g1_below_lod.below_detection_limit
    @test !g1_below_lod.abiotic_dominated
    @test g1_below_lod.biotic_mic_attribution_verdict == :Incomplete
    g1_missing = evaluate_cbm_g1_abiotic_compaction_gate(o2_excluded_anoxic = false)
    @test g1_missing.biotic_mic_attribution_verdict == :Incomplete
    @test_throws ArgumentError evaluate_cbm_g1_abiotic_compaction_gate(
        assume_higher_biomass_or_lower_density_means_higher_corrosion = true
    )

    g2_no_rheo = evaluate_cbm_g2_screening_vs_bridging_gate(
        proposed_mode = :polymer_bridging,
        require_empirical_rheology = false,
        measured_rheology_or_retention_shift = true
    )
    @test g2_no_rheo.verdict == :Incomplete
    g2_missing_meas = evaluate_cbm_g2_screening_vs_bridging_gate(
        proposed_mode = :polymer_bridging,
        require_empirical_rheology = true,
        measured_rheology_or_retention_shift = false
    )
    @test g2_missing_meas.verdict == :Incomplete
    g2_pass = evaluate_cbm_g2_screening_vs_bridging_gate(
        proposed_mode = :polymer_bridging,
        matched_monovalent_control = true,
        require_empirical_rheology = true,
        measured_rheology_or_retention_shift = true
    )
    @test g2_pass.verdict == :Accepted
    @test_throws ArgumentError evaluate_cbm_g2_screening_vs_bridging_gate(infer_bridging_from_zeta_potential_only = true)
    @test_throws ArgumentError evaluate_cbm_g2_screening_vs_bridging_gate(equate_prescribed_crosslink_fn_with_validation = true)

    g3_fe_only = evaluate_cbm_g3_diet_vs_epili_gate(magnetite_present = true)
    @test g3_fe_only.diet_causal_verdict == :Incomplete
    g3_no_pili = evaluate_cbm_g3_diet_vs_epili_gate(
        magnetite_present = true,
        e_pili_present_and_verified = false,
        direct_flux_perturbation_verified = true
    )
    @test g3_no_pili.diet_causal_verdict == :Rejected
    g3_verified = evaluate_cbm_g3_diet_vs_epili_gate(
        magnetite_present = true,
        e_pili_present_and_verified = true,
        soluble_shuttles_excluded = true,
        direct_flux_perturbation_verified = true
    )
    @test g3_verified.diet_causal_verdict == :Accepted
    @test_throws ArgumentError evaluate_cbm_g3_diet_vs_epili_gate(claim_diet_from_fe3o4_presence_alone = true)

    # Identifiability & Bead-Calibrated FCM Operator
    id_ok = evaluate_cbm_identifiability_and_fcm_operator()
    @test id_ok.identifiable
    @test id_ok.fisher_rank == 4
    @test id_ok.fcm_observation_operator.verdict == :Accepted
    @test !id_ok.fcm_observation_operator.resolves_clay_pore_spatial_position
    @test !id_ok.fcm_observation_operator.measures_eps_crosslink_density
    id_no_aw = evaluate_cbm_identifiability_and_fcm_operator(measure_aw_and_d001 = false)
    id_no_fe = evaluate_cbm_identifiability_and_fcm_operator(measure_fe_speciation = false)
    id_no_avs = evaluate_cbm_identifiability_and_fcm_operator(measure_sulfide_and_avs = false)
    @test !id_no_aw.identifiable && id_no_aw.fisher_rank < 4
    @test !id_no_fe.identifiable && id_no_fe.fisher_rank < 4
    @test !id_no_avs.identifiable && id_no_avs.fisher_rank < 4
    id_low = evaluate_cbm_identifiability_and_fcm_operator(
        measure_aw_and_d001 = false,
        measure_sulfide_and_avs = false
    )
    @test !id_low.identifiable
    @test id_low.fisher_rank < 4
    @test_throws ArgumentError evaluate_cbm_identifiability_and_fcm_operator(use_uncalibrated_raw_fcm_events = true)
    @test_throws ArgumentError evaluate_cbm_identifiability_and_fcm_operator(claim_spatial_position_from_fcm = true)

    # SYN6 & CBM Bundle JSON exports
    syn_b = export_syntrophic_eda_bundle_json()
    cbm_b = export_cbm_eda_bundle_json()
    @test isfile(syn_b.bundle_path)
    @test isfile(cbm_b.bundle_path)
    syn_txt = read(syn_b.bundle_path, String)
    cbm_txt = read(cbm_b.bundle_path, String)
    @test !occursin('\u2014', syn_txt) && !occursin('\u2013', syn_txt)
    @test !occursin('\u2014', cbm_txt) && !occursin('\u2013', cbm_txt)
    @test VALID_M7_STAGE_HEADERS == ("M7.1", "M7.2", "M7.3", "M7.4", "M7.5", "M7.6", "M7.7", "M7.8")
end
