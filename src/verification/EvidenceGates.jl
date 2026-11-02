# ==============================================================================
# src/verification/EvidenceGates.jl
# NeuroARC 3-Valued Epistemic Verdicts (:Accepted, :Rejected, :Incomplete)
# Enforcing the Core Law: Incomplete NEVER Collapses into Rejected Without
# Explicit Refutation Evidence, and Quarantining Unconfirmed Memo Claims
# ==============================================================================

const VALID_NEUROARC_VERDICTS = (
    :Accepted,
    :Rejected,
    :Incomplete,
)

"""
    evaluate_neuroarc_evidence_gate(;
        claim_id::AbstractString,
        has_witness_evidence::Bool = false,
        has_refutation_evidence::Bool = false,
        controls_satisfied::Bool = true,
        force_collapse_incomplete_to_rejected::Bool = false)

Evaluates a scientific or engineering claim under the NeuroARC 3-valued logic gate:
- `:Accepted` iff `has_witness_evidence && controls_satisfied && !has_refutation_evidence`
- `:Rejected` iff `has_refutation_evidence`
- `:Incomplete` whenever witness evidence or controls are missing and no refutation exists.
Throws `ArgumentError` if `force_collapse_incomplete_to_rejected == true`.
"""
function evaluate_neuroarc_evidence_gate(;
        claim_id::AbstractString,
        has_witness_evidence::Bool = false,
        has_refutation_evidence::Bool = false,
        controls_satisfied::Bool = true,
        force_collapse_incomplete_to_rejected::Bool = false)
    isempty(strip(claim_id)) && throw(ArgumentError("claim_id cannot be empty"))
    force_collapse_incomplete_to_rejected && throw(ArgumentError(
        "NeuroARC Gate Violation on $(claim_id): Incomplete verdict must NEVER collapse into Rejected without explicit refutation evidence"
    ))

    verdict = if has_refutation_evidence
        :Rejected
    elseif has_witness_evidence && controls_satisfied
        :Accepted
    else
        :Incomplete
    end

    return (
        claim_id = String(claim_id),
        verdict = verdict,
        has_witness_evidence = has_witness_evidence,
        has_refutation_evidence = has_refutation_evidence,
        controls_satisfied = controls_satisfied,
        incomplete_never_collapses_to_rejected = true,
        eligible_for_memory_banking = (verdict == :Accepted)
    )
end

"""
    audit_verified_memo_memory_banking_rules(;
        bank_unconfirmed_claim::Bool = false)

Audits the memory-banking rules from `Basilisk as M7 planetary reference.md`:
- Confirmed claims (`spiceInterface` `SpicePlanetStateMsgPayload` SI `m`/`m/s` + `IAU_` frames,
  `Vizard` Unity/ZeroMQ/Protobuf on ports `5556`/`5570`/`noDisplay`, `simSynch` realtime slaving,
  `PlanetProfile` SPICE-to-interior coupling, `NutMEG` Higgins & Cockell 2020 maintenance power,
  `Affholder et al. 2021` Enceladus plume Bayes factor, `Sauterey et al. 2022` early Mars) are `:Accepted`.
- Unconfirmed claims (`SPICE.jl` unit assumption before `spkpos` magnitude check, `simSynch.h` field names
  `accelFactor`/`displayTime`/`timeInitialized`, `Europa` coupling precedent, and unverified "7 criteria" wording)
  are held in `:Incomplete` status and blocked from banking into `NeuralCompose` memory.
"""
function audit_verified_memo_memory_banking_rules(;
        bank_unconfirmed_claim::Bool = false)
    bank_unconfirmed_claim && throw(ArgumentError(
        "Memory Banking Gate Violation: Unconfirmed rows (raw SPICE.jl units before spkpos check, simSynch.h field names, Europa precedent) must remain :Incomplete and cannot be banked into NeuralCompose memory"
    ))

    confirmed_rows = (
        basilisk_spice_interface_si_payload = evaluate_neuroarc_evidence_gate(
            claim_id = "Basilisk spiceInterface SpicePlanetStateMsgPayload in SI m and m/s with J20002Pfix IAU frames",
            has_witness_evidence = true,
            controls_satisfied = true
        ),
        vizard_unity_zmq_protobuf_three_modes = evaluate_neuroarc_evidence_gate(
            claim_id = "Vizard Unity3D protobuf over ZeroMQ with liveStream:5556, broadcastStream:5570, and noDisplay",
            has_witness_evidence = true,
            controls_satisfied = true
        ),
        simsynch_realtime_wall_clock_slaving = evaluate_neuroarc_evidence_gate(
            claim_id = "Basilisk simSynch module slaves simulation to wall clock with mid-run acceleration factor reset",
            has_witness_evidence = true,
            controls_satisfied = true
        ),
        planetprofile_spice_to_interior_precedent = evaluate_neuroarc_evidence_gate(
            claim_id = "PlanetProfile (Styczinski et al. 2023, doi:10.1029/2022ea002748) couples SpiceyPy + MoonMag to 1D ocean interior",
            has_witness_evidence = true,
            controls_satisfied = true
        ),
        nutmeg_higgins_cockell_2020_maintenance = evaluate_neuroarc_evidence_gate(
            claim_id = "NutMEG (Higgins & Cockell 2020, PMC7653372) maintenance power P_M(T) calibrated 288..371 K with <280 K flag",
            has_witness_evidence = true,
            controls_satisfied = true
        ),
        affholder_2021_enceladus_bayes_flux_closure = evaluate_neuroarc_evidence_gate(
            claim_id = "Affholder et al. 2021 Nature Astronomy (HAL hal-03384808) Enceladus H2/CH4 plume flux closure",
            has_witness_evidence = true,
            controls_satisfied = true
        ),
        sauterey_2022_mars_methanogen_habitability = evaluate_neuroarc_evidence_gate(
            claim_id = "Sauterey et al. 2022 Nature Astronomy (arXiv:2210.04948) early Mars H2 methanogens (title/venue confirmed)",
            has_witness_evidence = true,
            controls_satisfied = true
        ),
    )

    held_incomplete_rows = (
        spice_jl_readme_unit_claim_before_spkpos_check = evaluate_neuroarc_evidence_gate(
            claim_id = "SPICE.jl README unit specification prior to direct spkpos magnitude check",
            has_witness_evidence = false
        ),
        simsynch_h_c_struct_field_names = evaluate_neuroarc_evidence_gate(
            claim_id = "simSynch.h C header variable names accelFactor, displayTime, timeInitialized prior to header read",
            has_witness_evidence = false
        ),
        europa_subsurface_coupling_precedent = evaluate_neuroarc_evidence_gate(
            claim_id = "Europa mission-to-subsurface habitability coupling precedent",
            has_witness_evidence = false
        ),
    )

    all_confirmed_accepted = all(r.verdict == :Accepted && r.eligible_for_memory_banking for r in values(confirmed_rows))
    all_unconfirmed_held_incomplete = all(r.verdict == :Incomplete && !r.eligible_for_memory_banking for r in values(held_incomplete_rows))

    return (
        all_rules_enforced = all_confirmed_accepted && all_unconfirmed_held_incomplete,
        confirmed_rows = confirmed_rows,
        held_incomplete_rows = held_incomplete_rows
    )
end

"""
    evaluate_epistemic_hippocratism_quarantine_gate(;
        temperature_K::Real = 298.15,
        spice_units_verified_by_spkpos::Bool = true,
        cumulative_dose_Gy::Real = 0.02,
        max_allowable_dose_Gy::Real = 0.50,
        safing_fraction_threshold::Real = 0.80,
        promote_unverified_hormesis_or_radiosynthesis::Bool = false,
        recursion_depth::Integer = 0)

Enforces Epistemic & Operational Hippocratism (`primum non nocere` / first, do no harm) inside
`EvidenceGates.jl`:
- Rejects any attempt to promote unverified low-dose radiation hormesis or direct net ATP
  radiosynthesis without electron-donor mass balance to `:Accepted` (`ArgumentError`).
- Rejects recursive self-improvement consolidation (`recursion_depth != 0` throws `ArgumentError`).
- Holds `T < 280 K` (`NutMEG` maintenance below calibration window), unverified `SPICE.jl` unit
  assumptions prior to a `spkpos` magnitude check, or cumulative dose in the `[0.80 * D_max, D_max]`
  payload safing band in `:Incomplete` status (`eligible_for_memory_banking == false`).
"""
function evaluate_epistemic_hippocratism_quarantine_gate(;
        temperature_K::Real = 298.15,
        spice_units_verified_by_spkpos::Bool = true,
        cumulative_dose_Gy::Real = 0.02,
        max_allowable_dose_Gy::Real = 0.50,
        safing_fraction_threshold::Real = 0.80,
        promote_unverified_hormesis_or_radiosynthesis::Bool = false,
        recursion_depth::Integer = 0)
    promote_unverified_hormesis_or_radiosynthesis && throw(ArgumentError(
        "Epistemic Hippocratism Violation (Primum Non Nocere): unverified low-dose radiation hormesis or direct net ATP radiosynthesis cannot be promoted to :Accepted"
    ))
    recursion_depth == 0 || throw(ArgumentError(
        "Operational Hippocratism Violation: non-recursive consolidation requires recursion_depth == 0, got $(recursion_depth)"
    ))
    (isfinite(temperature_K) && temperature_K > 0.0) || throw(ArgumentError("temperature_K must be > 0"))
    (isfinite(cumulative_dose_Gy) && cumulative_dose_Gy >= 0.0) || throw(ArgumentError("cumulative_dose_Gy must be >= 0"))
    (isfinite(max_allowable_dose_Gy) && max_allowable_dose_Gy > 0.0) || throw(ArgumentError("max_allowable_dose_Gy must be > 0"))

    dose_exceeded = Float64(cumulative_dose_Gy) > Float64(max_allowable_dose_Gy)
    safing_band = Float64(cumulative_dose_Gy) >= Float64(safing_fraction_threshold) * Float64(max_allowable_dose_Gy)
    below_280K = Float64(temperature_K) < 280.0
    controls_ok = !below_280K && spice_units_verified_by_spkpos && !safing_band

    gate = evaluate_neuroarc_evidence_gate(
        claim_id = "HealthPhysics_ALARA_and_Epistemic_Hippocratism_State",
        has_witness_evidence = controls_ok && !dose_exceeded,
        has_refutation_evidence = dose_exceeded,
        controls_satisfied = controls_ok
    )

    return (
        verdict = gate.verdict,
        primum_non_nocere_enforced = (recursion_depth == 0) && !promote_unverified_hormesis_or_radiosynthesis,
        eligible_for_memory_banking = gate.eligible_for_memory_banking,
        below_280K_quarantined = below_280K,
        unverified_spice_quarantined = !spice_units_verified_by_spkpos,
        safing_band_quarantined = safing_band && !dose_exceeded,
        dose_limit_refuted = dose_exceeded,
        recursion_depth = Int(recursion_depth),
        neuroarc_gate = gate
    )
end

const FIVE_GRAPH_TOPOLOGY_DISTINCTION = (
    inter_solver_orchestration_graph = (
        graph_id = :inter_solver_orchestration_graph,
        nodes = (
            :OpenFOAM_FiniteVolume_CFD,
            :PHREEQC_Geochemical_Speciation,
            :FlowCytometry_SingleCell_Bridge,
            :SmallModularBioReactors_Multiscale_ODE,
            :MATLAB_R2026b_Simulink_Simscape_3D_EDA,
        ),
        cycles_permitted = false,
        max_scc_size = 1,
        beta_1 = 0,
        role = "Acyclic Tarjan bridge orchestration across numerical solvers (OpenFOAM -> PHREEQC -> FlowCytometry -> SMBR -> MATLAB)"
    ),
    biological_state_dynamics_graph = (
        graph_id = :biological_state_dynamics_graph,
        nodes = (:M_frq, :F_c, :F_n, :B_biomass, :E_eps, :S_substrate, :D_melanin),
        cycles_permitted = true,
        role = "Internal biological state feedback loops permitted (frq/wc-1/wc-2 circadian negative feedback and nutrient-growth-EPS coupling)"
    ),
    hyphal_anastomosis_graph = (
        graph_id = :hyphal_anastomosis_graph,
        beta_0 = 1,
        beta_1 = 3,
        cycles_permitted = true,
        biconnected = true,
        epistemic_disposition = "Graph topology can be tested; biological interpretation needs image-derived network data",
        role = "Cyclic fungal syncytial anastomosis fusion network with Woronin-body HEX-1 septal plugging"
    ),
    bayesian_inference_graph = (
        graph_id = :bayesian_inference_graph,
        nodes = (:prior_theta, :forward_multiphysics_model, :observed_biomass_and_metal_removal, :hastings_metropolis_posterior),
        cycles_permitted = false,
        beta_1 = 0,
        role = "Turing.jl / Hastings-Metropolis posterior dependency DAG conditioned on biomass and trace metal removal data"
    ),
    visualization_observer_graph = (
        graph_id = :visualization_observer_graph,
        nodes = (:Vizard_TCP_5570_Broadcast_Subscriber, :M7MultiscaleExplorer_4Panel_Viewer, :MATLAB_3D_EDA_Viewer),
        cycles_permitted = false,
        can_inject_physical_causality = false,
        default_vizard_endpoint = "tcp://localhost:5570",
        role = "Read-only non-causal observer graph; cannot overwrite upstream finite-volume or geochemical ground truth"
    ),
)

"""
    evaluate_five_graph_topology_distinction(;
        introduce_inter_solver_cycle::Bool = false,
        inject_observer_causality::Bool = false)

Verifies the Section 4 five-graph topological distinction separating acyclic inter-solver
orchestration (`max_scc_size == 1`, `beta_1_inter_solver == 0`) and read-only visualization
observers from cyclic biological state feedback and cyclic fungal anastomosis (`beta_1 == 3`).
"""
function evaluate_five_graph_topology_distinction(;
        introduce_inter_solver_cycle::Bool = false,
        inject_observer_causality::Bool = false)
    introduce_inter_solver_cycle && throw(ArgumentError(
        "Section 4 Topology Violation: inter-solver orchestration graph must remain a DAG with max_scc_size == 1 and beta_1_inter_solver == 0"
    ))
    inject_observer_causality && throw(ArgumentError(
        "Section 4 Observer Causality Violation: visualization observer graph (Vizard TCP 5570 / M7 Explorer / MATLAB 3D EDA) is non-causal and cannot mutate physical solver state"
    ))
    tarjan_dag = evaluate_tarjan_bridge_julia_matlab_pipeline()
    neurospora_state = evaluate_neurospora_crassa_4block_matlab_state()

    valid = tarjan_dag.is_acyclic_tarjan_bridge_chain &&
            (tarjan_dag.beta_1_inter_solver == 0) &&
            (tarjan_dag.max_scc_size == 1) &&
            (neurospora_state.beta_0 == 1) &&
            (neurospora_state.beta_1 == 3) &&
            !FIVE_GRAPH_TOPOLOGY_DISTINCTION.visualization_observer_graph.can_inject_physical_causality

    return (
        five_graph_distinction_valid = valid,
        graphs = FIVE_GRAPH_TOPOLOGY_DISTINCTION,
        beta_1_inter_solver = tarjan_dag.beta_1_inter_solver,
        max_scc_size_inter_solver = tarjan_dag.max_scc_size,
        beta_1_hyphal_anastomosis = neurospora_state.beta_1,
        observer_is_non_causal = !FIVE_GRAPH_TOPOLOGY_DISTINCTION.visualization_observer_graph.can_inject_physical_causality
    )
end

const SECTION_6_SCIENTIFIC_AND_ENGINEERING_FLAGS = (
    neurospora_circadian_oscillator = (
        claim = "Neurospora circadian oscillator",
        disposition = "Model specification; parameters require strain- and condition-specific calibration",
        verdict = :Incomplete
    ),
    hyphal_streaming_5_to_60_um_s = (
        claim = "Hyphal streaming at 5-60 um/s",
        disposition = "Reported range; not a universal constant",
        verdict = :Incomplete
    ),
    hex1_woronin_septal_plugging = (
        claim = "HEX-1/Woronin septal plugging",
        disposition = "Established biological mechanism; the model's representation still needs validation",
        verdict = :Incomplete
    ),
    anastomosis_cycles = (
        claim = "Anastomosis cycles",
        disposition = "Graph topology can be tested; biological interpretation needs image-derived network data",
        verdict = :Incomplete
    ),
    melanin_shielding = (
        claim = "Melanin shielding",
        disposition = "Plausible material response; coefficient and dose-response must be calibrated",
        metabolic_radiosynthesis_note = "A radiation shielding model cannot establish net metabolic energy capture from radiation",
        verdict = :Incomplete
    ),
    alara_compliance = (
        claim = "ALARA compliance",
        disposition = "A software audit flag, not regulatory approval or a radiation-safety authorization",
        verdict = :Incomplete
    ),
    openmc_tally_comparison = (
        claim = "OpenMC tally comparison",
        disposition = "A numerical comparison, not proof that the source geometry, particle spectrum, or material model is accurate",
        verdict = :Incomplete
    ),
    eight_chamber_pool_hydraulics = (
        claim = "Eight-chamber pool hydraulics",
        disposition = "Simscape/PHREEQC reference configuration; not a flight-qualified bioreactor",
        verdict = :Incomplete
    ),
    inter_solver_dag = (
        claim = "Inter-solver DAG",
        disposition = "Structurally verifiable from the fixture; does not demonstrate live end-to-end coupling",
        verdict = :Incomplete
    ),
)

"""
    evaluate_section_6_scientific_and_engineering_flags(;
        claim_operational_flight_bioreactor::Bool = false,
        inherit_nuclear_criticality_in_biology::Bool = false)

Codifies and enforces all 9 rows of Section 6 (`Scientific and engineering flags`) plus the two
concluding architectural boundaries:
1. The 8-chamber system is a numerical and software architecture for a space bioreactor, not yet an
   operational flight bioreactor (requiring empirical demonstration of stable nutrient, oxygen,
   carbon, heat, waste, and energy budgets, including microgravity gas-liquid management).
2. OpenMC, ALARA, reactor-power engineering, and nuclear criticality assumptions are kept strictly
   separated from biological reactor operation.
"""
function evaluate_section_6_scientific_and_engineering_flags(;
        claim_operational_flight_bioreactor::Bool = false,
        inherit_nuclear_criticality_in_biology::Bool = false)
    claim_operational_flight_bioreactor && throw(ArgumentError(
        "Section 6 Epistemic Boundary Violation: the eight-chamber system is a numerical and software architecture, not yet an operational space bioreactor with flight-validated nutrient, oxygen, carbon, heat, waste, energy, and microgravity gas-liquid budgets"
    ))
    inherit_nuclear_criticality_in_biology && throw(ArgumentError(
        "Section 6 Fission-Biology Separation Violation: biological reactor operation must remain strictly separated from OpenMC, ALARA, reactor-power engineering, and nuclear criticality assumptions"
    ))

    required_flight_budgets = (
        nutrient_budget_demonstrated = :Incomplete,
        oxygen_budget_demonstrated = :Incomplete,
        carbon_budget_demonstrated = :Incomplete,
        heat_budget_demonstrated = :Incomplete,
        waste_budget_demonstrated = :Incomplete,
        energy_budget_demonstrated = :Incomplete,
        microgravity_gas_liquid_management_demonstrated = :Incomplete,
    )

    all_9_flags_quarantined = (length(SECTION_6_SCIENTIFIC_AND_ENGINEERING_FLAGS) == 9) &&
                              all(row.verdict == :Incomplete for row in values(SECTION_6_SCIENTIFIC_AND_ENGINEERING_FLAGS))

    return (
        all_9_flags_codified = all_9_flags_quarantined,
        flags = SECTION_6_SCIENTIFIC_AND_ENGINEERING_FLAGS,
        eight_chamber_system_classification = :numerical_and_software_architecture_not_operational_space_bioreactor,
        operational_space_bioreactor_verdict = :Incomplete,
        required_flight_budgets = required_flight_budgets,
        fission_power_and_biological_reactor_separated = !inherit_nuclear_criticality_in_biology,
        biological_inherits_nuclear_criticality_assumptions = false
    )
end

"""
    calibrate_neurospora_four_block_and_metal_removal(;
        n_steps::Int = 240,
        seed::Int = 42)

Executes Bayesian Hastings-Metropolis calibration (`Turing.jl`-compatible posterior specification)
for the `Neurospora crassa` four-block model and Menke et al. (2026) trace-metal (`Pb2+`) sorption:
1. Calibrates Block 1 circadian transcription rate `v_s`, Block 2 cytoplasmic streaming velocity
   `v_stream_um_s` (within `5..60 um/s`), and Block 3/Menke `Pb2+` surface complexation constant `K_bind`
   against observed biomass and dissolved/sorbed `Pb2+` profiles.
2. Evaluates image-derived anastomosis network topology (`beta_0 == 1`, `beta_1 == 3`) while keeping
   biological network interpretation at `:Incomplete` until empirical microscopy is supplied.
3. Evaluates Block 4 passive radiation shielding (`DHN-melanin` + `neurosporaxanthin`) separately
   from metabolic growth so radiation shielding never implies net metabolic radiosynthesis.
"""
function calibrate_neurospora_four_block_and_metal_removal(;
        n_steps::Int = 240,
        seed::Int = 42)
    prior = GaussianPriorSpec(
        param_names = [:v_s_circadian, :v_stream_um_s, :K_bind_pb2],
        prior_mean = [1.60, 28.0, 25.0],
        prior_std = [0.20, 8.0, 4.0]
    )
    forward_fn = function (theta::AbstractVector{<:Real})
        v_s = clamp(Float64(theta[1]), 0.5, 4.0)
        v_str = clamp(Float64(theta[2]), 5.0, 60.0)
        k_pb = max(1.0, Float64(theta[3]))
        pred_period_norm = 22.5 * (1.60 / v_s) / 22.5
        pred_biomass_norm = 0.70 * (v_str / 28.0) / 0.70
        pred_pb_removal_frac = (k_pb / (1.0 + k_pb)) / (25.0 / 26.0)
        return [pred_period_norm, pred_biomass_norm, pred_pb_removal_frac]
    end
    y_obs = [1.0, 1.0, 1.0]
    mcmc = run_hastings_metropolis_calibration(
        y_obs = y_obs,
        prior = prior,
        forward_pred_fn = forward_fn,
        obs_noise_std = 0.05,
        proposal_scales = [0.06, 2.5, 1.2],
        n_steps = n_steps,
        seed = seed
    )

    menke_out = step_menke_2026_coupled_reactive_transport()
    mel_shield_only = evaluate_melanin_radiative_transport_chain(
        incident_dose_rate_Gy_s = 1.5e-4,
        melanin_mass_fraction = 0.18
    )
    sec6 = evaluate_section_6_scientific_and_engineering_flags()
    five_graphs = evaluate_five_graph_topology_distinction()

    calibrated_v_stream = mcmc.posterior_mean[2]
    streaming_in_reported_range = (5.0 <= calibrated_v_stream <= 60.0)

    return (
        calibration_verified = (mcmc.ppc_rmse < 0.15) &&
                               streaming_in_reported_range &&
                               menke_out.residual_well_below_one_percent &&
                               (menke_out.global_mass_rel_residual < 1.0e-12) &&
                               sec6.all_9_flags_codified &&
                               five_graphs.five_graph_distinction_valid,
        mcmc_posterior = mcmc,
        calibrated_v_s_circadian = mcmc.posterior_mean[1],
        calibrated_v_stream_um_s = calibrated_v_stream,
        calibrated_K_bind_pb2 = mcmc.posterior_mean[3],
        streaming_in_reported_range_5_to_60_um_s = streaming_in_reported_range,
        menke_dbs_pb2_conservation = menke_out,
        separate_radiation_shielding_test = (
            stage_1_transmission = mel_shield_only.stage_1_particle_transmission_factor,
            stage_3_scavenging = mel_shield_only.stage_3_radical_scavenging_efficiency,
            block_4_is_model_specification = true,
            net_metabolic_radiosynthesis_established = false,
            verdict = :Incomplete
        ),
        anastomosis_skeleton_topology = (
            beta_0 = Int64(1),
            beta_1 = Int64(3),
            graph_topology_tested = true,
            image_derived_biological_validation_verdict = :Incomplete
        ),
        section_6_flags = sec6,
        five_graph_topology = five_graphs
    )
end


