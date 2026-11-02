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

