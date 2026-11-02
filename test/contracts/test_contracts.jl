@testset "Contracts: UnitSystem, EvidenceClass, SpatialFrame, ReactorState" begin
    # 1. UnitSystem
    q1 = SMBRUnitQuantity("pore_velocity", 1.5e-4, :m_s)
    q2 = SMBRUnitQuantity("pressure", 101325.0, :Pa)
    hom = verify_unit_homogeneity([q1, q2])
    @test hom.homogeneous
    @test hom.quantity_count == 2
    @test_throws ArgumentError SMBRUnitQuantity("", 1.0, :m)
    @test_throws ArgumentError SMBRUnitQuantity("bad_inf", Inf, :m)
    @test_throws ArgumentError SMBRUnitQuantity("bad_unit", 1.0, :furlongs)
    @test_throws ArgumentError verify_unit_homogeneity(SMBRUnitQuantity[])

    si_pos = convert_naif_km_to_si_m([1.495978707e8, 0.0, -250.0])
    @test isapprox(si_pos[1], 1.495978707e11; rtol = 1.0e-12)
    @test isapprox(si_pos[3], -2.5e5; rtol = 1.0e-12)
    @test_throws ArgumentError convert_naif_km_to_si_m([1.0, 2.0])
    @test_throws ArgumentError convert_naif_km_to_si_m([NaN, 0.0, 0.0])

    c_molar = convert_molality_to_molarity_si(0.002, 1000.0; tds_kg_kg = 0.0)
    @test isapprox(c_molar, 2.0; atol = 1.0e-12)
    @test_throws ArgumentError convert_molality_to_molarity_si(-0.1, 1000.0)
    @test_throws ArgumentError convert_molality_to_molarity_si(0.01, 0.0)
    @test_throws ArgumentError convert_molality_to_molarity_si(0.01, 1000.0; tds_kg_kg = 1.0)

    ph_res = compute_activity_averaged_ph([7.35, 6.85, 7.35, 6.85])
    @test isapprox(ph_res.ph_activity_avg, 7.03169896; atol = 1.0e-6)
    @test isapprox(ph_res.ph_arithmetic_invalid, 7.10; atol = 1.0e-9)
    @test ph_res.arithmetic_vs_activity_bias > 0.05
    @test_throws ArgumentError compute_activity_averaged_ph(Float64[])
    @test_throws ArgumentError compute_activity_averaged_ph([15.0])

    # 2. EvidenceClass & Role Separation
    r_solv = ArchitecturalRoleContract(component_name = "OpenFOAMSolver", assigned_role = :physical_solver, mutates_physical_state = true)
    @test r_solv.assigned_role == :physical_solver
    @test_throws ArgumentError ArchitecturalRoleContract(
        component_name = "ConflatedSolver",
        assigned_role = :physical_solver,
        mutates_physical_state = true,
        performs_statistical_inference = true
    )
    @test_throws ArgumentError ArchitecturalRoleContract(
        component_name = "BadVisualizer",
        assigned_role = :visualization,
        mutates_physical_state = true
    )

    img_man = verify_downloads_five_images_manifest()
    @test img_man.manifest_count == 5
    @test img_man.all_hashes_match
    @test img_man.fission_vs_bioreactor_strictly_separated
    @test img_man.disk_verified_count == 5

    # 3. SpatialFrame & Reversible Coordinate Transform
    tf = SMBRCoordinateTransform(
        source_frame = :reactor_local_um,
        target_frame = :spacecraft_payload_m,
        scale_to_meters = 1.0e-6,
        translation_m = [0.12, -0.05, 0.20]
    )
    pt_um = [250.0, -100.0, 50.0]
    pt_m = transform_point_forward(tf, pt_um)
    pt_rec = transform_point_inverse(tf, pt_m)
    @test isapprox(pt_rec, pt_um; atol = 1.0e-9)
    @test_throws ArgumentError SMBRCoordinateTransform(rotation_matrix = [-1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])

    frame_ok = SMBRSpatialFrame()
    fv_ok = validate_spatial_frame_for_fv_solver(frame_ok; require_conservative_mesh = true)
    @test fv_ok.valid && fv_ok.fv_eligible
    frame_splat = SMBRSpatialFrame(mesh_kind = :gaussian_splat, is_watertight_manifold = false, closed_face_fluxes = false)
    @test_throws ArgumentError validate_spatial_frame_for_fv_solver(frame_splat; require_conservative_mesh = true)

    # 4. ReactorState, 3-Fraction Ledger, Global Conservation, Clay Transitions
    ledger = ThreeFractionAnalyticalLedger(
        element_symbol = :Pb,
        aqueous_dissolved_mol = 1.0e-6,
        mineral_adsorbed_mol = 2.0e-6,
        cellular_polymer_mol = 0.5e-6
    )
    @test isapprox(total_element_moles(ledger), 3.5e-6; atol = 1.0e-15)
    @test ledger.icp_ms_instrument == :ICP_MS
    @test ledger.xrd_xrf_instrument == :XRD_XRF
    @test ledger.flow_cytometry_instrument == :Flow_Cytometry

    rs = ReactorState()
    @test rs.cell_count == 4
    @test rs.clay_system == :system_c_illite_smectite_eps_biofilm
    @test_throws ArgumentError ReactorState(phi = fill(1.2, 4))

    cons_ok = evaluate_global_mass_conservation_rate(
        initial_total_mol = 1.0e-4,
        final_total_mol = 1.05e-4,
        dt_s = 10.0,
        J_in_mol_s = 1.0e-6,
        J_out_mol_s = 0.5e-6,
        R_nuclear_mol_s = 0.0
    )
    @test cons_ok.is_conserved && cons_ok.well_below_one_percent
    cons_bad = evaluate_global_mass_conservation_rate(
        initial_total_mol = 1.0e-4,
        final_total_mol = 1.5e-4,
        dt_s = 10.0,
        J_in_mol_s = 0.0,
        J_out_mol_s = 0.0
    )
    @test !cons_bad.is_conserved

    ca = evaluate_clay_system_transition(:system_a_sudoite_illite)
    cb = evaluate_clay_system_transition(:system_b_bentonite_montmorillonite)
    cc = evaluate_clay_system_transition(:system_c_illite_smectite_eps_biofilm)
    @test cb.cec_meq_100g > cc.cec_meq_100g > ca.cec_meq_100g
    @test_throws ArgumentError evaluate_clay_system_transition(:unknown_clay)
end
