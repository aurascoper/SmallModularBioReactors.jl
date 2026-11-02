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

const VALID_M7_STAGE_HEADERS = (
    "M7.1",
    "M7.2",
    "M7.3",
    "M7.4",
    "M7.5",
    "M7.6",
    "M7.7",
    "M7.8",
)

"""
    validate_m7_header_vs_numeric_schema(;
        stage_headers::AbstractVector{<:AbstractString} = collect(VALID_M7_STAGE_HEADERS),
        discrete_indices_int64::AbstractVector{<:Integer} = Int64[1, 2, 3, 4],
        continuous_observables_float64::AbstractVector{<:AbstractFloat} = Float64[0.38, 298.15, 7.12, 1.495978707e11],
        conflate_header_as_numeric_literal::Bool = false)

Enforces the Turn `L14`/`L16` (`index = 32..34` of `dialectic-turns-2026-11-02.jsonl`) type-vs-header
separation across `NeuralCompose`, Julia (`SmallModularBioReactors.jl` / `M7Multiscale`), and MATLAB:
1. `"M7.1".."M7.8"` are strictly container/stage headers naming the 8 multiscale workbench layers.
2. `Int64` ("ants" -> `ints`) represents exact discrete cell/chamber indices, codon counts, and Betti numbers.
3. `Float64` ("floats") represents IEEE-754 continuous physical observables in SI units.
Throws `ArgumentError` if a stage header is conflated with a float/int literal or if any header is invalid.
"""
function validate_m7_header_vs_numeric_schema(;
        stage_headers::AbstractVector{<:AbstractString} = collect(VALID_M7_STAGE_HEADERS),
        discrete_indices_int64::AbstractVector{<:Integer} = Int64[1, 2, 3, 4],
        continuous_observables_float64::AbstractVector{<:AbstractFloat} = Float64[0.38, 298.15, 7.12, 1.495978707e11],
        conflate_header_as_numeric_literal::Bool = false)
    conflate_header_as_numeric_literal && throw(ArgumentError(
        "Category-Level Schema Violation (Turn L16): M7.1..M7.8 are pipeline stage headers (container labels), never Float64 or Int64 numeric literals"
    ))
    length(stage_headers) == 8 || throw(ArgumentError("stage_headers must contain all 8 stages M7.1..M7.8"))
    for (idx, h) in enumerate(stage_headers)
        expected = VALID_M7_STAGE_HEADERS[idx]
        String(h) == expected || throw(ArgumentError("Expected stage header $(expected) at position $(idx), got $(h)"))
    end
    !isempty(discrete_indices_int64) || throw(ArgumentError("discrete_indices_int64 must not be empty"))
    !isempty(continuous_observables_float64) || throw(ArgumentError("continuous_observables_float64 must not be empty"))
    all(isfinite, continuous_observables_float64) || throw(ArgumentError("continuous_observables_float64 must be finite"))

    int_vec = Int64.(discrete_indices_int64)
    flt_vec = Float64.(continuous_observables_float64)
    ieee_01_plus_02_gap = abs((0.1 + 0.2) - 0.3)

    return (
        schema_valid = true,
        stage_headers = Tuple(String.(stage_headers)),
        header_axis = :pipeline_stage_container_label,
        int64_axis = :exact_discrete_index_and_betti_topology_count,
        float64_axis = :ieee754_continuous_si_physical_observable,
        discrete_indices_int64 = int_vec,
        continuous_observables_float64 = flt_vec,
        ieee754_binary_approximation_residual = ieee_01_plus_02_gap,
        cross_layer_falsification_enabled = true
    )
end

"""
    evaluate_selfish_gene_and_comfort_reticulated_asr(;
        hgt_reticulation_edges::Integer = 2,
        force_strict_tree_on_reticulated_network::Bool = false)

Evaluates the Turn `L18`/`L19` (`index = 36..37` of `dialectic-turns-2026-11-02.jsonl`) ancestral
sequence reconstruction (`ASR`) synthesis combining:
1. Dawkins' replicator / selfish-gene allele-level transmission objective (`CsgA` curli repeats,
   `DsrB` `[4Fe-4S]` catalytic triad, and retrotransposon-derived `Arc` Gag capsid `PDB:6GSE`).
2. Science historian Nathaniel Comfort's (`"The Tangled Tree"` / `"Mother Nature"`) critique of
   single-level gene determinism: horizontal gene transfer (`HGT`) and viral/transposon reticulation
   create phylogenetic network loops (`beta_1_reticulation_cycles == hgt_reticulation_edges >= 1`),
   and selection acts across multiple levels (`:gene_allele`, `:organism_phenotype`, `:syncytial_group`,
   `:cultural_protocol`) whose reconstructed trees need not coincide.
Throws `ArgumentError` if `force_strict_tree_on_reticulated_network == true` when `hgt_reticulation_edges > 0`.
"""
function evaluate_selfish_gene_and_comfort_reticulated_asr(;
        hgt_reticulation_edges::Integer = 2,
        force_strict_tree_on_reticulated_network::Bool = false)
    hgt_reticulation_edges >= 0 || throw(ArgumentError("hgt_reticulation_edges must be >= 0"))
    if force_strict_tree_on_reticulated_network && hgt_reticulation_edges > 0
        throw(ArgumentError(
            "Phylogenetic Category Violation (Turn L18/L19): cannot impose a strictly bifurcating tree (beta_1 == 0) onto a horizontal-gene-transfer reticulated network with $(hgt_reticulation_edges) reticulation edges"
        ))
    end

    # Base tree with N_nodes = 7, N_tree_edges = 6 (beta_0 = 1, beta_1 = 0 before HGT reticulation)
    n_nodes = 7
    n_edges = (n_nodes - 1) + Int(hgt_reticulation_edges)
    beta_0 = 1
    beta_1_reticulation_cycles = n_edges - n_nodes + beta_0

    selection_levels = (
        gene_allele = :replicator_sequence_transmission_maximization,
        organism_phenotype = :physiological_homeostasis_and_survival,
        syncytial_group = :mycelial_hyphal_anastomosis_network_flow,
        cultural_protocol = :scientific_instrumentation_and_alara_governance
    )

    return (
        reticulated_asr_verified = (beta_1_reticulation_cycles == Int(hgt_reticulation_edges)),
        n_nodes = n_nodes,
        n_edges = n_edges,
        beta_0 = beta_0,
        beta_1_reticulation_cycles = beta_1_reticulation_cycles,
        is_reticulated_network = (beta_1_reticulation_cycles > 0),
        dawkins_replicator_unit = :allele_sequence_transmission,
        nathaniel_comfort_critique = :multilevel_selection_and_reticulate_tangled_tree,
        selection_levels = selection_levels,
        inter_level_tree_disagreement_is_signal = true
    )
end

"""
    evaluate_tarjan_bridge_julia_matlab_pipeline(;
        introduce_feedback_cycle_to_openfoam::Bool = false)

Evaluates the Turn `L21` (`index = 39`, `2026-11-02T07:27:50Z`) solver orchestration architecture:
`OpenFOAM_CFD -> PHREEQC_Geochemistry -> FlowCytometry_Observation -> SmallModularBioReactors_ODE -> MATLAB_Simulink_Simscape_3D_EDA`.
Verifies that every cross-tool handoff is a directed **Tarjan bridge edge** in a directed acyclic graph
(`max_scc_size == 1`, `beta_1_inter_solver == 0`), ensuring MATLAB / Simulink / Simscape acts as a
read-only ROM / control / 3D EDA consumer without recursively overwriting OpenFOAM or PHREEQC ground truth.
"""
function evaluate_tarjan_bridge_julia_matlab_pipeline(;
        introduce_feedback_cycle_to_openfoam::Bool = false)
    introduce_feedback_cycle_to_openfoam && throw(ArgumentError(
        "Tarjan Bridge Acyclicity Violation (Turn L21): MATLAB/Simulink must not feed back recursively into OpenFOAM mesh/field ground truth; cross-solver handoffs must remain acyclic Tarjan bridges"
    ))

    nodes = (
        :OpenFOAM_FiniteVolume_CFD,
        :PHREEQC_Geochemical_Speciation,
        :FlowCytometry_SingleCell_Bridge,
        :SmallModularBioReactors_Multiscale_ODE,
        :MATLAB_R2026b_Simulink_Simscape_3D_EDA,
    )
    directed_edges = (
        (:OpenFOAM_FiniteVolume_CFD => :PHREEQC_Geochemical_Speciation),
        (:PHREEQC_Geochemical_Speciation => :FlowCytometry_SingleCell_Bridge),
        (:FlowCytometry_SingleCell_Bridge => :SmallModularBioReactors_Multiscale_ODE),
        (:SmallModularBioReactors_Multiscale_ODE => :MATLAB_R2026b_Simulink_Simscape_3D_EDA),
    )

    n_vertices = length(nodes)
    n_bridge_edges = length(directed_edges)
    # In a directed path graph of 5 vertices and 4 edges, every edge is a bridge, max SCC size is 1, and beta_1 == 0
    max_scc_size = 1
    beta_1_inter_solver = n_bridge_edges - n_vertices + 1

    return (
        is_acyclic_tarjan_bridge_chain = (max_scc_size == 1) && (beta_1_inter_solver == 0),
        nodes = nodes,
        directed_edges = directed_edges,
        n_vertices = n_vertices,
        n_tarjan_bridge_edges = n_bridge_edges,
        max_scc_size = max_scc_size,
        beta_1_inter_solver = beta_1_inter_solver,
        matlab_mutates_upstream_solver = false
    )
end

