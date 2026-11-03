@testset "Biology & Molecular: Biofilm CA, NutMEG 7 Criteria, Radiotropism, Melanin, Sobolev" begin
    # 1. Conservative CA Biofilm Redistribution (arXiv:2610.01565)
    ca_step = step_conservative_ca_biofilm_redistribution(
        B_attached_kg_m3 = [1.35, 0.50, 0.40, 0.30],
        B_suspended_kg_m3 = [0.10, 0.08, 0.05, 0.04],
        phi = fill(0.38, 4),
        B_max_cell_kg_m3 = 1.20
    )
    @test ca_step.biomass_conserved
    @test ca_step.B_attached_next_kg_m3[1] <= 1.20 + 1.0e-12

    ca_multi = step_conservative_ca_biofilm_redistribution(
        B_attached_kg_m3 = [1.18, 1.35, 0.40, 0.30],
        B_suspended_kg_m3 = zeros(4),
        phi = fill(0.38, 4),
        B_max_cell_kg_m3 = 1.20,
        attachment_rate_s_inv = 0.0,
        detachment_rate_s_inv = 0.0
    )
    @test ca_multi.biomass_conserved
    @test all(x -> x <= 1.20 + 1.0e-12, ca_multi.B_attached_next_kg_m3)

    # 2. Section 5 Traceable 7 Budget-Closure Criteria (NutMEG, Affholder 2021, Sauterey 2022, PlanetProfile)
    nm_280 = evaluate_nutmeg_methanogen_budget_closure(temperature_K = 280.0, target_world = :ENCELADUS)
    nm_330 = evaluate_nutmeg_methanogen_budget_closure(temperature_K = 330.0, target_world = :ENCELADUS)
    @test !nm_280.out_of_calibrated_maintenance_range_below_280K
    @test nm_330.within_calibrated_288_to_371K
    @test isapprox(nm_280.final_biomass_cells_L, 1.0e9; rtol = 1.0e-6)
    @test isapprox(nm_330.final_biomass_cells_L, 1.0e6; rtol = 1.0e-6)
    @test isapprox(nm_280.total_ch4_uM, 62.0; rtol = 1.0e-6)
    @test isapprox(nm_330.total_ch4_uM, 3.0; rtol = 1.0e-6)
    @test isapprox(nm_280.peak_growth_rate_h_inv, 0.007; rtol = 1.0e-6)
    @test isapprox(nm_330.peak_growth_rate_h_inv, 0.14; rtol = 1.0e-6)

    # Below 280 K out-of-range flag
    nm_cold = evaluate_nutmeg_methanogen_budget_closure(temperature_K = 273.15, target_world = :ENCELADUS)
    @test nm_cold.out_of_calibrated_maintenance_range_below_280K
    @test nm_cold.world_verdict == :Incomplete

    # 330 K with 90% of optimal CO2 yields near-zero growth
    nm_co2_limited = evaluate_nutmeg_methanogen_budget_closure(
        temperature_K = 330.0,
        co2_fraction_of_optimal = 0.90,
        dissolved_phosphorus_mol_L = 1.0e-10
    )
    @test nm_co2_limited.peak_growth_rate_h_inv < 0.001

    # Degenerate fixed-point ATP yield must throw ArgumentError
    @test_throws ArgumentError evaluate_nutmeg_methanogen_budget_closure(
        atp_yield_interval_mol_per_mol_co2 = (1.0, 1.0)
    )

    # Sauterey 2022 assigned to Mars (not Enceladus) and Europa held in :Incomplete
    @test nm_330.precedent_attribution.mars_habitability_and_cooling_feedback.target_planet == :MARS
    @test !nm_330.precedent_attribution.mars_habitability_and_cooling_feedback.applies_to_enceladus
    @test nm_330.precedent_attribution.europa_coupling_precedent.verdict == :Incomplete
    @test !nm_330.planetprofile_coupling.habitat_volume_is_free_parameter

    # 3. RNAReplicatorToyAdapter Quarantine
    toy0 = RNAReplicatorToyState()
    @test toy0.epistemic_tag == :synthetic_2d_lattice_toy_model
    @test !toy0.is_physical_3d_pore_solver
    toy1 = step_rna_replicator_toy(toy0)
    @test toy1.epistemic_tag == :synthetic_2d_lattice_toy_model
    @test_throws ArgumentError RNAReplicatorToyState(promote_to_3d_continuum_pore_solver = true)

    # 4. Radiotropism, Gravitropism, MelaninResponse
    rt_ok = evaluate_radiotropism_flux_and_controls(
        controls = RadiotropismControlBattery(true, true, true, true, false)
    )
    @test rt_ok.verdict == :Accepted
    rt_inc = evaluate_radiotropism_flux_and_controls(
        controls = RadiotropismControlBattery(true, false, true, true, false)
    )
    @test rt_inc.verdict == :Incomplete
    rt_rej = evaluate_radiotropism_flux_and_controls(
        controls = RadiotropismControlBattery(true, false, true, true, true)
    )
    @test rt_rej.verdict == :Rejected

    gf = compute_gravitropic_flux()
    @test gf.flux_magnitude_kg_m2_s > 0.0

    mel = evaluate_melanin_radiative_transport_chain(melanin_mass_fraction = 0.20)
    @test mel.stage_3_radical_scavenging_efficiency > 0.40
    @test mel.quarantined_hypotheses.melanin_mediated_magnetoreception == :unestablished_hypothesis
    @test_throws ArgumentError evaluate_melanin_radiative_transport_chain(enable_unestablished_magnetoreception = true)
    @test_throws ArgumentError evaluate_melanin_radiative_transport_chain(enable_unestablished_net_atp_radiosynthesis = true)

    nc_preset = NeurosporaCrassaCircadianHyphalPreset()
    @test nc_preset.organism_symbol == :Neurospora_crassa
    @test nc_preset.clock_genes == (:frq, :wc_1, :wc_2)
    @test nc_preset.free_running_period_hr == 22.5
    nc_dd = evaluate_neurospora_crassa_syncytial_circadian_preset(nc_preset; elapsed_time_hr = 2.0, optical_entrainment_active = false)
    nc_ld = evaluate_neurospora_crassa_syncytial_circadian_preset(nc_preset; elapsed_time_hr = 2.0, optical_entrainment_active = true)
    nc_zero_rad = evaluate_neurospora_crassa_syncytial_circadian_preset(nc_preset; elapsed_time_hr = 2.0, incident_dose_rate_Gy_s = 0.0)
    @test nc_dd.active_period_hr == 22.5
    @test nc_ld.active_period_hr == 24.0
    @test nc_dd.syncytial_signal_propagation_intact
    @test length(nc_dd.chamber_hyphal_extension_um_hr) == 4
    @test nc_zero_rad.mean_hyphal_extension_um_hr > 100.0
    @test 5.0 < nc_dd.mean_hyphal_extension_um_hr < nc_zero_rad.mean_hyphal_extension_um_hr
    @test_throws ArgumentError NeurosporaCrassaCircadianHyphalPreset(free_running_period_hr = 5.0)
    @test_throws ArgumentError evaluate_neurospora_crassa_syncytial_circadian_preset(nc_preset; enable_unestablished_magnetoreception = true)

    # 5. Molecular: SobolevAdapter, ProteinPeptideRegistry, StructureProvenance
    anc = SobolevMolecularAnchor(sobolev_relaxation_step = 200)
    raw_sm = encode_smbr_smacro1_bytes(anc.sobolev_relaxation_step, anc.bead_coordinates_nm)
    dec_sm = decode_smbr_smacro1_bytes(raw_sm)
    @test dec_sm.sobolev_relaxation_step == 200
    @test isapprox(dec_sm.bead_coordinates_nm, anc.bead_coordinates_nm; atol = 1.0e-12)
    @test_throws ArgumentError decode_smbr_smacro1_bytes(UInt8[0x00, 0x01])

    csga = lookup_protein_peptide_entry(:CsgA)
    dsrb = lookup_protein_peptide_entry(:DsrB)
    arc_cap = lookup_protein_peptide_entry(:Arc)
    @test csga.core_length_aa == 12 && csga.repeat_length_aa == 18
    @test occursin("n_X = 5", dsrb.motif_signature)
    @test arc_cap.uniprot_or_pdb_id == "PDB:6GSE" && arc_cap.stoichiometric_binding_sites_per_monomer == 240
    @test_throws ArgumentError lookup_protein_peptide_entry(:NonExistentProtein)

    sp = StructureAndBindingProvenance()
    @test sp.provenance_separated
    @test_throws ArgumentError StructureAndBindingProvenance(claim_kd_from_static_pdb_coordinates_only = true)
    @test_throws ArgumentError StructureAndBindingProvenance(binding_constant_K_bind_m3_mol = 10.0, binding_assay_method = nothing)

    rna_led = verify_rna_biofilm_toy_mass_ledger()
    @test rna_led.mass_ledger_closed
    @test rna_led.residual <= 1.0e-8
    @test rna_led.diffusion_clip == 0.0
    @test_throws ArgumentError verify_rna_biofilm_toy_mass_ledger(deaths = 2, death_return = 25.0)

    bf_gate = evaluate_biofilms_cpm_calibration_refusal_gate()
    @test bf_gate.refusal_contracts_preserved
    @test bf_gate.clock == :numerical_mcs
    @test isnan(bf_gate.seconds_per_mcs)
    @test bf_gate.lattice_pitch_m === nothing
    @test bf_gate.radiodialysis_coupled_gate == "RADIODIALYSIS: BLOCKED"
    @test bf_gate.evidence_status == :Incomplete
    @test_throws ArgumentError evaluate_biofilms_cpm_calibration_refusal_gate(attempt_radiodialysis_on_occupancy_mean = true)
end

