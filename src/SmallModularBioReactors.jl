module SmallModularBioReactors

using Dates
using LinearAlgebra
using Printf
using Random
using SHA
import SPICE

@inline function mean(x::AbstractArray; dims = nothing)
    if dims === nothing
        return sum(x) / length(x)
    elseif dims == 2
        return reshape([sum(@view x[i, :]) / size(x, 2) for i in 1:size(x, 1)], :, 1)
    elseif dims == 1
        return reshape([sum(@view x[:, j]) / size(x, 1) for j in 1:size(x, 2)], 1, :)
    else
        throw(ArgumentError("Unsupported dims=$(dims) in mean"))
    end
end

@inline function var(x::AbstractVector)
    n = length(x)
    n <= 1 && return 0.0
    m = mean(x)
    return sum((Float64(v) - m)^2 for v in x) / (n - 1)
end

@inline function var(x::AbstractMatrix; dims = nothing)
    if dims === nothing
        return var(vec(x))
    elseif dims == 2
        return reshape([var(@view x[i, :]) for i in 1:size(x, 1)], :, 1)
    elseif dims == 1
        return reshape([var(@view x[:, j]) for j in 1:size(x, 2)], 1, :)
    else
        throw(ArgumentError("Unsupported dims=$(dims) in var"))
    end
end

# 1. Contracts
include(joinpath(@__DIR__, "contracts", "UnitSystem.jl"))
include(joinpath(@__DIR__, "contracts", "EvidenceClass.jl"))
include(joinpath(@__DIR__, "contracts", "SpatialFrame.jl"))
include(joinpath(@__DIR__, "contracts", "ReactorState.jl"))

# 2. Environment
include(joinpath(@__DIR__, "environment", "MissionClock.jl"))
include(joinpath(@__DIR__, "environment", "SPICEAdapter.jl"))
include(joinpath(@__DIR__, "environment", "GravityField.jl"))
include(joinpath(@__DIR__, "environment", "RadiationField.jl"))

# 3. Spacecraft
include(joinpath(@__DIR__, "spacecraft", "BasiliskAdapter.jl"))
include(joinpath(@__DIR__, "spacecraft", "PowerBudget.jl"))
include(joinpath(@__DIR__, "spacecraft", "ThermalControl.jl"))
include(joinpath(@__DIR__, "spacecraft", "PayloadController.jl"))

# 4. Transport
include(joinpath(@__DIR__, "transport", "OpenFOAMAdapter.jl"))
include(joinpath(@__DIR__, "transport", "PHREEQCAdapter.jl"))
include(joinpath(@__DIR__, "transport", "PolymerTransport.jl"))
include(joinpath(@__DIR__, "transport", "MineralExchange.jl"))

# 5. Biology
include(joinpath(@__DIR__, "biology", "BiofilmDynamics.jl"))
include(joinpath(@__DIR__, "biology", "RNAReplicatorToyAdapter.jl"))
include(joinpath(@__DIR__, "biology", "Radiotropism.jl"))
include(joinpath(@__DIR__, "biology", "Gravitropism.jl"))
include(joinpath(@__DIR__, "biology", "MelaninResponse.jl"))

# 6. Molecular
include(joinpath(@__DIR__, "molecular", "SobolevAdapter.jl"))
include(joinpath(@__DIR__, "molecular", "ProteinPeptideRegistry.jl"))
include(joinpath(@__DIR__, "molecular", "StructureProvenance.jl"))

# 7. Inference
include(joinpath(@__DIR__, "inference", "BayesianCalibration.jl"))
include(joinpath(@__DIR__, "inference", "HastingsMetropolis.jl"))
include(joinpath(@__DIR__, "inference", "DriftDetection.jl"))
include(joinpath(@__DIR__, "inference", "ParetoDesign.jl"))

# 8. Visualization
include(joinpath(@__DIR__, "visualization", "M7Adapter.jl"))
include(joinpath(@__DIR__, "visualization", "VizardAdapter.jl"))

# 9. Verification
include(joinpath(@__DIR__, "verification", "Conservation.jl"))
include(joinpath(@__DIR__, "verification", "EvidenceGates.jl"))
include(joinpath(@__DIR__, "inference", "EvidenceGates.jl"))
include(joinpath(@__DIR__, "verification", "Uncertainty.jl"))
include(joinpath(@__DIR__, "verification", "AcceptanceTests.jl"))

export VALID_SI_BASE_UNITS,
       NAIF_KM_TO_SI_M,
       SMBRUnitQuantity,
       convert_naif_km_to_si_m,
       convert_molality_to_molarity_si,
       compute_activity_averaged_ph,
       verify_unit_homogeneity,
       VALID_ARCHITECTURAL_ROLES,
       VALID_EPISTEMIC_TIERS,
       ArchitecturalRoleContract,
       DOWNLOADS_FIVE_IMAGES_MANIFEST,
       verify_downloads_five_images_manifest,
       VALID_SMBR_MESH_KINDS,
       SMBRCoordinateTransform,
       transform_point_forward,
       transform_point_inverse,
       SMBRSpatialFrame,
       validate_spatial_frame_for_fv_solver,
       VALID_CLAY_SYSTEMS,
       VALID_ANALYTICAL_INSTRUMENTS,
       ThreeFractionAnalyticalLedger,
       total_element_moles,
       ReactorState,
       VALID_REACTOR_POOL_IDS,
       MultiPoolReactorBiofilmState,
       evaluate_global_mass_conservation_rate,
       evaluate_clay_system_transition,
       MissionClock,
       advance_mission_clock,
       SMBR_SPICE_LOCK,
       STANDARD_SPICE_KERNELS,
       VALID_SPICE_TARGETS,
       BASILISK_SPICE_INTERFACE_OUTPUT_MSGS,
       BASILISK_SPICE_INTERFACE_METHODS,
       SpicePlanetStateMsgPayload,
       convert_spice_state_to_basilisk_si,
       SPICEKernelSet,
       spice_furnsh,
       spice_utc2et,
       spice_pxform,
       spice_spkpos,
       execute_live_spice_jl_known_body_check,
       CANONICAL_SURFACE_GRAVITY_M_S2,
       compute_effective_gravity,
       SpectralRadiationFieldSpec,
       compute_dose_rate_profile,
       HealthPhysicsALARAContract,
       evaluate_health_physics_alara_and_hippocratism,
       evaluate_openmc_cross_section_substrate_swap_gate,
       BasiliskSpacecraftState,
       step_basilisk_spacecraft_state,
       ExternalPowerSupplySpec,
       BioreactorPowerDemand,
       evaluate_power_budget,
       step_reactor_thermal_balance,
       PayloadSafetyControllerSpec,
       evaluate_alara_payload_controller_safing,
       compute_kozeny_carman_permeability,
       step_darcy_brinkman_stokes_transport,
       DEFAULT_OPENFOAM_POLYMESH_DIR,
       parse_openfoam_polymesh_topology_and_roundtrip,
       PHREEQCCellSpeciation,
       evaluate_phreeqc_cell_array,
       MICROPILLAR_STREAMER_TIMEPOINTS_HR,
       simulate_micropillar_streamer_bridging,
       solve_polymer_ligand_metal_complexation,
       step_menke_2026_coupled_reactive_transport,
       TRACEABLE_SEVEN_BUDGET_CRITERIA,
       step_conservative_ca_biofilm_redistribution,
       evaluate_nutmeg_methanogen_budget_closure,
       coordinate_matlab_multipool_biofilm_session_seam,
       RNAReplicatorToyState,
       step_rna_replicator_toy,
       RadiotropismControlBattery,
       evaluate_radiotropism_flux_and_controls,
       compute_gravitropic_flux,
       QUARANTINED_MELANIN_HYPOTHESES,
       evaluate_melanin_radiative_transport_chain,
       NeurosporaCrassaCircadianHyphalPreset,
       evaluate_neurospora_crassa_syncytial_circadian_preset,
       evaluate_melanin_neurospora_alara_shielding_factor,
       evaluate_neurospora_crassa_4block_matlab_state,
       SMACRO1_MAGIC_BYTES,
       SobolevMolecularAnchor,
       encode_smbr_smacro1_bytes,
       decode_smbr_smacro1_bytes,
       ProteinPeptideEntry,
       PROTEIN_PEPTIDE_REGISTRY,
       lookup_protein_peptide_entry,
       StructureAndBindingProvenance,
       GaussianPriorSpec,
       evaluate_log_posterior,
       compute_metropolis_hastings_alpha,
       run_hastings_metropolis_calibration,
       PosteriorDriftRecord,
       compute_posterior_drift,
       SIX_PARETO_OBJECTIVE_NAMES,
       ParetoDesignCandidate,
       pareto_dominates,
       extract_pareto_front,
       CRITERIAL_STATUS_CODES,
       build_m7_workbench_scene,
       VALID_M7_STAGE_HEADERS,
       validate_m7_header_vs_numeric_schema,
       evaluate_selfish_gene_and_comfort_reticulated_asr,
       evaluate_tarjan_bridge_julia_matlab_pipeline,
       evaluate_matlab_neuralcompose_collar_and_multipool_biofilm_seam,
       VIZARD_ZMQ_HANDSHAKE_STRINGS,
       VIZARD_LIVE_MODES,
       VizardTransportConfig,
       build_vizard_vizmessage_packet,
       audit_reactor_conservation_invariants,
       VALID_NEUROARC_VERDICTS,
       evaluate_neuroarc_evidence_gate,
       audit_verified_memo_memory_banking_rules,
       evaluate_epistemic_hippocratism_quarantine_gate,
       FIVE_GRAPH_TOPOLOGY_DISTINCTION,
       evaluate_five_graph_topology_distinction,
       SECTION_6_SCIENTIFIC_AND_ENGINEERING_FLAGS,
       evaluate_section_6_scientific_and_engineering_flags,
       calibrate_neurospora_four_block_and_metal_removal,
       decompose_epistemic_aleatoric_uncertainty,
       audit_equifinality_underdetermination,
       run_demonstration_1_clay_eps_metal_reactor,
       run_demonstration_2_radiation_gradient_phenotype,
       run_demonstration_3_spacecraft_mission_reactor,
       evaluate_smbr_milestones_s0_to_s6

end # module SmallModularBioReactors
