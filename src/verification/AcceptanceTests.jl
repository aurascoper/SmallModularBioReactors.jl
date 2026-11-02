# ==============================================================================
# src/verification/AcceptanceTests.jl
# End-to-End Runners for Demonstrations 1, 2, 3 and Milestones S0..S6
# ==============================================================================

"""
    run_demonstration_1_clay_eps_metal_reactor()

Executes Demonstration 1: Clay-EPS-Metal Microfluidic Reactor (`System A/B/C` clay transition,
`Menke et al. 2026` conservative pore-water + adsorbed surface complexation `P-L + M^{2+} <=> P-L-M`,
microfluidic micropillar streamer bridging `t = 5..13 hr` anchored to `rOkmBjSK...jpg`, and
3-fraction `ICP-MS` + `XRD/XRF` + `Flow Cytometry` ledger).
"""
function run_demonstration_1_clay_eps_metal_reactor()
    clay_a = evaluate_clay_system_transition(:system_a_sudoite_illite)
    clay_b = evaluate_clay_system_transition(:system_b_bentonite_montmorillonite)
    clay_c = evaluate_clay_system_transition(:system_c_illite_smectite_eps_biofilm; eps_fraction = 0.22)

    menke_out = step_menke_2026_coupled_reactive_transport()
    streamer_out = simulate_micropillar_streamer_bridging()
    ledger = ThreeFractionAnalyticalLedger(
        element_symbol = :Pb,
        aqueous_dissolved_mol = sum(menke_out.C_aq_mol_m3) * 0.38e-9,
        mineral_adsorbed_mol = sum(menke_out.M_sorb_mol_m3_bulk) * 1.0e-9,
        cellular_polymer_mol = 1.2e-10
    )

    passed = menke_out.residual_well_below_one_percent &&
             (menke_out.global_mass_rel_residual < 1.0e-12) &&
             streamer_out.streamer_bridging_established &&
             (total_element_moles(ledger) > 0.0) &&
             (clay_b.cec_meq_100g > clay_c.cec_meq_100g > clay_a.cec_meq_100g)

    return (
        demonstration_id = :Demonstration_1_Clay_EPS_Metal_Microfluidic_Reactor,
        passed = passed,
        clay_transitions = (system_a = clay_a, system_b = clay_b, system_c = clay_c),
        menke_reactive_transport = menke_out,
        micropillar_streamers = streamer_out,
        three_fraction_ledger = ledger
    )
end

"""
    run_demonstration_2_radiation_gradient_phenotype()

Executes Demonstration 2: Radiation-Gradient Phenotype & Melanin Response (`RadiationField` spectral
integral `D_dot(x,t)` and gradient `grad(D_dot)`, `MelaninResponse` 4-stage transport chain comparing
melanized vs non-melanized control colonies anchored to `4I3OYyXt...jpg`, and `Radiotropism` flux with
mandatory 4-control artifact battery).
"""
function run_demonstration_2_radiation_gradient_phenotype()
    rad_spec = SpectralRadiationFieldSpec()
    depths = [0.0, 0.002, 0.004, 0.006]
    rad_profile = compute_dose_rate_profile(rad_spec, depths; melanin_mass_fraction = 0.18)

    melanized = evaluate_melanin_radiative_transport_chain(
        incident_dose_rate_Gy_s = rad_profile.dose_rate_Gy_s[1],
        melanin_mass_fraction = 0.20
    )
    unmelanized_control = evaluate_melanin_radiative_transport_chain(
        incident_dose_rate_Gy_s = rad_profile.dose_rate_Gy_s[1],
        melanin_mass_fraction = 0.0
    )

    tropism_controlled = evaluate_radiotropism_flux_and_controls(
        dose_rate_Gy_s = rad_profile.dose_rate_Gy_s,
        grad_dose_rate_Gy_s_m = rad_profile.grad_dose_rate_Gy_s_m,
        controls = RadiotropismControlBattery(true, true, true, true, false)
    )
    tropism_uncontrolled = evaluate_radiotropism_flux_and_controls(
        dose_rate_Gy_s = rad_profile.dose_rate_Gy_s,
        grad_dose_rate_Gy_s_m = rad_profile.grad_dose_rate_Gy_s_m,
        controls = RadiotropismControlBattery(false, true, true, true, false)
    )

    neurospora_syncytium = evaluate_neurospora_crassa_syncytial_circadian_preset(
        NeurosporaCrassaCircadianHyphalPreset();
        elapsed_time_hr = 6.0,
        optical_entrainment_active = true,
        incident_dose_rate_Gy_s = rad_profile.dose_rate_Gy_s[1]
    )

    passed = rad_profile.spice_geometry_decoupled &&
             (melanized.stage_3_unscavenged_radical_rate_mol_m3_s < unmelanized_control.stage_3_unscavenged_radical_rate_mol_m3_s) &&
             (melanized.stage_4_net_radial_growth_um_hr > unmelanized_control.stage_4_net_radial_growth_um_hr) &&
             (tropism_controlled.verdict == :Accepted) &&
             (tropism_uncontrolled.verdict == :Incomplete) &&
             neurospora_syncytium.syncytial_signal_propagation_intact &&
             (neurospora_syncytium.mean_hyphal_extension_um_hr > 0.0)

    return (
        demonstration_id = :Demonstration_2_Radiation_Gradient_Phenotype,
        passed = passed,
        radiation_profile = rad_profile,
        melanized_colony = melanized,
        unmelanized_control_colony = unmelanized_control,
        radiotropism_controlled = tropism_controlled,
        radiotropism_uncontrolled = tropism_uncontrolled,
        neurospora_crassa_syncytium = neurospora_syncytium
    )
end

"""
    run_demonstration_3_spacecraft_mission_reactor()

Executes Demonstration 3: Spacecraft Mission-Linked Reactor (`SPICE.jl` + `Basilisk` `SpicePlanetStateMsgPayload`
via single boundary converter `convert_spice_state_to_basilisk_si`, ESA `BIM` microgravity carousel
centrifuge `ytJAGCW0...jpg`, `Vizard` `broadcastStream` on `tcp://localhost:5570` + `M7Adapter`,
Section 9 external Kilopower fission supply vs internal Ambr 250 biological reactor separation, and
Section 5 Traceable 7 Budget-Closure Criteria).
"""
function run_demonstration_3_spacecraft_mission_reactor()
    kernels = spice_furnsh()
    et0 = spice_utc2et("2026-11-01T12:00:00Z"; kernels = kernels)
    sc0 = BasiliskSpacecraftState(
        scenario_name = "scenarioMonteCarloSpice_ESA_BIM",
        monte_carlo_run_index = 1,
        mission_elapsed_time_s = 0.0,
        target_body = :MARS,
        bim_carousel_rpm = 15.0
    )
    step_out = step_basilisk_spacecraft_state(sc0, 60.0; et0_s = et0, carousel_radius_m = 0.14)

    clock = advance_mission_clock(
        MissionClock();
        d_mission_s = 60.0,
        d_flow_s = 5.0,
        d_kinetic_s = 5.0,
        d_sobolev_steps = 10,
        sync_event_id = "SYNC_BIM_CAROUSEL_STEP_1"
    )

    power_check = evaluate_power_budget(ExternalPowerSupplySpec(), BioreactorPowerDemand())
    therm_check = step_reactor_thermal_balance()
    grav_flux = compute_gravitropic_flux(g_eff_m_s2 = step_out.effective_gravity.g_eff_m_s2)

    tf = SMBRCoordinateTransform()
    vizard_pkt = build_vizard_vizmessage_packet(step_out.spacecraft_state, clock, tf)
    m7_scene = build_m7_workbench_scene(ReactorState(), SMBRSpatialFrame(), clock; verdict = :Accepted)

    nutmeg_300K = evaluate_nutmeg_methanogen_budget_closure(temperature_K = 300.0, target_world = :ENCELADUS)
    nutmeg_275K_cold = evaluate_nutmeg_methanogen_budget_closure(temperature_K = 275.0, target_world = :ENCELADUS)
    nutmeg_europa = evaluate_nutmeg_methanogen_budget_closure(temperature_K = 300.0, target_world = :EUROPA)

    passed = power_check.budget_satisfied &&
             power_check.fission_supply_external_only &&
             power_check.bioreactor_is_not_fission_simulator &&
             therm_check.within_biocompatible_window &&
             (grav_flux.flux_magnitude_kg_m2_s > 0.0) &&
             (vizard_pkt.active_mode == :broadcastStream) &&
             (vizard_pkt.transform_roundtrip_error_um <= 1.0e-9) &&
             (m7_scene.criterial_status_code == 1.0) &&
             !nutmeg_300K.out_of_calibrated_maintenance_range_below_280K &&
             nutmeg_275K_cold.out_of_calibrated_maintenance_range_below_280K &&
             (nutmeg_europa.world_verdict == :Incomplete)

    return (
        demonstration_id = :Demonstration_3_Spacecraft_Mission_Linked_Reactor,
        passed = passed,
        basilisk_step = step_out,
        power_budget = power_check,
        thermal_control = therm_check,
        gravitropic_flux = grav_flux,
        vizard_packet = vizard_pkt,
        m7_scene = m7_scene,
        nutmeg_300K = nutmeg_300K,
        nutmeg_275K_cold = nutmeg_275K_cold,
        nutmeg_europa = nutmeg_europa,
        evidence_image_sha256 = DOWNLOADS_FIVE_IMAGES_MANIFEST.demonstration_3_esa_bim_microgravity_carousel.sha256
    )
end

"""
    evaluate_smbr_milestones_s0_to_s6()

Runs the complete Milestone `S0..S6` acceptance verification suite for `SmallModularBioReactors.jl`.
"""
function evaluate_smbr_milestones_s0_to_s6()
    d1 = run_demonstration_1_clay_eps_metal_reactor()
    d2 = run_demonstration_2_radiation_gradient_phenotype()
    d3 = run_demonstration_3_spacecraft_mission_reactor()
    img_manifest = verify_downloads_five_images_manifest()
    memo_rules = audit_verified_memo_memory_banking_rules()
    mcmc = run_hastings_metropolis_calibration(n_steps = 120)
    drift = compute_posterior_drift(
        [0.12, 25.0, 1.5e-6],
        Matrix{Float64}(I, 3, 3) .* 0.01,
        mcmc.posterior_mean,
        mcmc.posterior_cov
    )
    pareto = extract_pareto_front([
        ParetoDesignCandidate(
            candidate_id = "design_A_optimal_balance",
            biomass_productivity = 0.85,
            metal_toxicity = 0.05,
            radiation_damage = 0.04,
            energy_consumption = 95.0,
            permeability_loss = 0.12,
            mass_and_volume = 4.2
        ),
        ParetoDesignCandidate(
            candidate_id = "design_B_dominated",
            biomass_productivity = 0.60,
            metal_toxicity = 0.15,
            radiation_damage = 0.10,
            energy_consumption = 120.0,
            permeability_loss = 0.30,
            mass_and_volume = 6.0
        )
    ])

    s0_passed = d3.power_budget.fission_supply_external_only &&
                d3.power_budget.bioreactor_is_not_fission_simulator
    s1_passed = d3.basilisk_step.spice_ephemeris.basilisk_payload.planet_frame == "IAU_MARS" &&
                d2.radiation_profile.spice_geometry_decoupled
    s2_passed = d1.passed
    s3_passed = d2.passed && (d3.nutmeg_europa.world_verdict == :Incomplete) &&
                d3.nutmeg_275K_cold.out_of_calibrated_maintenance_range_below_280K
    s4_passed = (lookup_protein_peptide_entry(:CsgA).core_length_aa == 12) &&
                (lookup_protein_peptide_entry(:DsrB).core_length_aa == 8) &&
                StructureAndBindingProvenance().provenance_separated
    s5_passed = mcmc.sampling_acceptance_is_not_engineering_verdict &&
                drift.separated_from_physical_state &&
                (pareto.pareto_front_count == 1)
    s6_passed = memo_rules.all_rules_enforced &&
                img_manifest.all_hashes_match &&
                img_manifest.fission_vs_bioreactor_strictly_separated &&
                d3.passed

    all_passed = s0_passed && s1_passed && s2_passed && s3_passed && s4_passed && s5_passed && s6_passed

    return (
        all_milestones_passed = all_passed,
        milestones = (
            S0_contracts_units_clocks_and_fission_boundary = s0_passed,
            S1_spice_basilisk_gravity_and_radiation_separation = s1_passed,
            S2_openfoam_phreeqc_menke2026_and_polymer_transport = s2_passed,
            S3_biofilm_nutmeg_radiotropism_and_melanin_quarantine = s3_passed,
            S4_sobolev_smacro1_and_structure_binding_provenance = s4_passed,
            S5_bayesian_mcmc_drift_and_pareto_inference = s5_passed,
            S6_neuroarc_gates_vizard_m7_and_downloads_manifest = s6_passed,
        ),
        demonstrations = (d1 = d1, d2 = d2, d3 = d3),
        downloads_manifest = img_manifest,
        memo_banking_audit = memo_rules
    )
end
