@testset "Environment, Spacecraft, and Porous Reactive Transport" begin
    # 1. MissionClock (Three-Clock Law)
    clk0 = MissionClock()
    clk1 = advance_mission_clock(clk0; d_mission_s = 120.0, d_flow_s = 10.0, d_kinetic_s = 5.0, d_sobolev_steps = 25)
    @test clk1.mission_elapsed_time_s == 120.0
    @test clk1.physical_flow_time_s == 10.0
    @test clk1.kinetic_biogeochem_time_s == 5.0
    @test clk1.sobolev_relaxation_step == 25
    @test !clk1.clocks_explicitly_coupled
    @test_throws ArgumentError advance_mission_clock(clk0; treat_sobolev_as_physical_time = true)
    @test_throws ArgumentError advance_mission_clock(clk0; assume_clocks_equal_without_event = true)
    @test_throws ArgumentError MissionClock(clocks_explicitly_coupled = true, sync_event_id = nothing)

    # 2. SPICEAdapter & Single Boundary Converter convert_spice_state_to_basilisk_si
    kset = spice_furnsh()
    @test kset.loaded
    @test kset.spice_jl_version == "0.2.5"
    et = spice_utc2et("2026-11-01T12:00:00Z"; kernels = kset)
    @test et > 8.0e8
    @test_throws ArgumentError spice_utc2et("invalid-date"; kernels = kset)

    # Unit test comparing convert_spice_state_to_basilisk_si against a Basilisk-generated SpicePlanetStateMsgPayload for the same epoch
    pos_spice_km = [2.279392e8, -1.500000e7, 4.200000e6]
    vel_spice_km_s = [1.580, 24.025, 0.510]
    dcm_iau = spice_pxform("J2000", "IAU_MARS", et)
    bsk_ref_payload = SpicePlanetStateMsgPayload(
        PositionVector = [2.279392e11, -1.500000e10, 4.200000e9],
        VelocityVector = [1580.0, 24025.0, 510.0],
        J2000Current = et,
        J20002Pfix = dcm_iau,
        PlanetName = "MARS"
    )
    conv_unchecked = convert_spice_state_to_basilisk_si(pos_spice_km, vel_spice_km_s, dcm_iau; j2000_current_s = et, planet_name = "MARS")
    @test conv_unchecked.spice_jl_unit_verdict == :Incomplete
    @test isapprox(conv_unchecked.payload.PositionVector, bsk_ref_payload.PositionVector; rtol = 1.0e-12)
    @test isapprox(conv_unchecked.payload.VelocityVector, bsk_ref_payload.VelocityVector; rtol = 1.0e-12)
    @test conv_unchecked.payload.J2000Current == bsk_ref_payload.J2000Current
    @test isapprox(conv_unchecked.payload.J20002Pfix, bsk_ref_payload.J20002Pfix; atol = 1.0e-12)
    @test conv_unchecked.payload.planet_frame == "IAU_MARS"

    conv_checked = convert_spice_state_to_basilisk_si(
        pos_spice_km,
        vel_spice_km_s,
        dcm_iau;
        j2000_current_s = et,
        planet_name = "MARS",
        expected_distance_m_bounds = (2.0e11, 2.6e11)
    )
    @test conv_checked.spice_jl_unit_verdict == :Accepted

    hst1 = spice_spkpos(:HST_LEO, et; monte_carlo_run_index = 1, verify_magnitude_bounds = true)
    hst2 = spice_spkpos(:HST_LEO, et; monte_carlo_run_index = 2, verify_magnitude_bounds = true)
    @test hst1.spice_jl_unit_verdict == :Accepted
    @test hst2.altitude_m > hst1.altitude_m
    @test !hst1.radiation_dose_transport_included
    @test_throws ArgumentError spice_spkpos(:MARS, et; request_radiation_dose_from_spice = true)

    # 3. GravityField & RadiationField
    g_mars = compute_effective_gravity(body = :MARS, omega_rad_s = [0.0, 0.0, 1.5], r_offset_m = [0.14, 0.0, 0.0])
    @test g_mars.g_eff_magnitude_m_s2 > 3.72
    @test_throws ArgumentError compute_effective_gravity(body = :MARS, g_planetary_m_s2 = [0.0, 0.0, -9.80665])

    rad_prof = compute_dose_rate_profile(SpectralRadiationFieldSpec(), [0.0, 0.002, 0.004, 0.006]; melanin_mass_fraction = 0.15)
    @test issorted(rad_prof.dose_rate_Gy_s; rev = true)
    @test all(g < 0.0 for g in rad_prof.grad_dose_rate_Gy_s_m)
    @test rad_prof.spice_geometry_decoupled

    # 4. Spacecraft: BasiliskAdapter, PowerBudget, ThermalControl
    sc_step = step_basilisk_spacecraft_state(BasiliskSpacecraftState(), 30.0; et0_s = et)
    @test sc_step.spacecraft_state.mission_elapsed_time_s == 30.0
    @test sc_step.carousel_omega_rad_s > 0.0

    pb = evaluate_power_budget(ExternalPowerSupplySpec(), BioreactorPowerDemand())
    @test pb.budget_satisfied
    @test pb.fission_supply_external_only
    @test pb.bioreactor_is_not_fission_simulator
    @test_throws ArgumentError ExternalPowerSupplySpec(is_external_to_bioreactor = false)
    @test_throws ArgumentError BioreactorPowerDemand(simulate_fission_criticality_inside_bioreactor = true)

    th = step_reactor_thermal_balance()
    @test th.within_biocompatible_window
    @test th.thermal_gradient_negligible_for_tropism

    # 5. Transport: OpenFOAM DBS, PHREEQC, PolymerTransport, MineralExchange (Menke et al. 2026)
    dbs = step_darcy_brinkman_stokes_transport(
        C_aq_mol_m3 = [0.30, 0.20, 0.10, 0.05],
        phi = fill(0.38, 4),
        closed_boundaries = true
    )
    @test dbs.mass_conserved_below_1pct
    @test dbs.rel_mass_residual < 1.0e-12
    @test_throws ArgumentError step_darcy_brinkman_stokes_transport(
        C_aq_mol_m3 = [0.30, 0.20],
        phi = [0.38, 0.38],
        u_darcy_m_s = 0.5,
        dx_m = 1.0e-4,
        dt_s = 10.0
    )

    ph_cells = [PHREEQCCellSpeciation(cell_index = i, pH = 7.0 + 0.1 * i) for i in 1:4]
    ph_aud = evaluate_phreeqc_cell_array(ph_cells)
    @test ph_aud.charge_balanced
    @test ph_aud.arithmetic_vs_activity_bias > 0.0

    strm = simulate_micropillar_streamer_bridging()
    @test strm.streamer_bridging_established
    @test last(strm.beta_1_anastomosis_loops) == 2

    eq_plm = solve_polymer_ligand_metal_complexation(
        total_metal_mol_m3_bulk = 0.25,
        total_ligand_mol_m3_bulk = 0.40,
        porosity = 0.38,
        K_bind_m3_mol = 30.0
    )
    @test eq_plm.mass_balance_rel_residual < 1.0e-14
    @test 0.0 < eq_plm.site_occupancy_fraction < 1.0

    menke = step_menke_2026_coupled_reactive_transport()
    @test menke.residual_well_below_one_percent
    @test menke.global_mass_rel_residual < 1.0e-12
end
