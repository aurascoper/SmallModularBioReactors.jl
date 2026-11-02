# ==============================================================================
# src/visualization/M7Adapter.jl
# Read-Only Visualization Bridge Exporting ReactorState + SMBRSpatialFrame +
# MissionClock to M7 4-Panel Workbench & SMACRO1 Molecular Anchors with
# Criterial Separation (+1.0 Accepted, 0.0 Incomplete, -1.0 Rejected)
# ==============================================================================

const CRITERIAL_STATUS_CODES = Dict{Symbol,Float64}(
    :Accepted => 1.0,
    :Incomplete => 0.0,
    :Rejected => -1.0,
)

"""
    build_m7_workbench_scene(state::ReactorState, frame::SMBRSpatialFrame, clock::MissionClock;
        verdict::Symbol = :Accepted,
        molecular_anchor::Union{Nothing,SobolevMolecularAnchor} = SobolevMolecularAnchor())

Constructs a read-only 4-panel M7 scientific workbench scene payload without mutating
physical solver state, preserving the Three-Clock Law and criterial separation codes.
"""
function build_m7_workbench_scene(
        state::ReactorState,
        frame::SMBRSpatialFrame,
        clock::MissionClock;
        verdict::Symbol = :Accepted,
        molecular_anchor::Union{Nothing,SobolevMolecularAnchor} = SobolevMolecularAnchor())
    haskey(CRITERIAL_STATUS_CODES, verdict) ||
        throw(ArgumentError("verdict must be one of :Accepted, :Incomplete, :Rejected"))
    state.cell_count == frame.cell_count ||
        throw(ArgumentError("ReactorState cell_count ($(state.cell_count)) must match SMBRSpatialFrame ($(frame.cell_count))"))

    fv_check = validate_spatial_frame_for_fv_solver(frame; require_conservative_mesh = false)
    status_code = CRITERIAL_STATUS_CODES[verdict]

    panel_a = (
        panel_id = :panel_a_3d_geometry_flow,
        mesh_kind = frame.mesh_kind,
        cell_count = frame.cell_count,
        mean_porosity = mean(state.phi),
        mean_permeability_m2 = mean(state.k_m2),
        beta_0 = frame.beta_0,
        beta_1 = frame.beta_1
    )
    panel_b = (
        panel_id = :panel_b_geochemistry,
        clay_system = state.clay_system,
        aqueous_species = collect(keys(state.C_aq_mol_m3)),
        mineral_surface_species = collect(keys(state.M_min_mol_m3))
    )
    panel_c = (
        panel_id = :panel_c_radiolytic_thermal,
        mean_temperature_K = mean(state.T_K),
        radiolytic_species = collect(keys(state.N_iso_mol_m3))
    )
    panel_d = (
        panel_id = :panel_d_biological_molecular,
        mean_biomass_kg_m3 = mean(state.B_kg_m3),
        mean_eps_kg_m3 = mean(state.E_kg_m3),
        mean_peptide_mol_m3 = mean(state.P_mol_m3),
        sobolev_anchor = molecular_anchor
    )

    return (
        role_contract = ArchitecturalRoleContract(
            component_name = "M7Adapter",
            assigned_role = :visualization,
            renders_scene_only = true
        ),
        verdict = verdict,
        criterial_status_code = status_code,
        clocks = clock,
        frame_sha256 = fv_check.frame_sha256,
        panels = (
            panel_a = panel_a,
            panel_b = panel_b,
            panel_c = panel_c,
            panel_d = panel_d
        )
    )
end
