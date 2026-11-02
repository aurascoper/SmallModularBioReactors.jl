# ==============================================================================
# src/spacecraft/PayloadController.jl
# Automated ALARA Time (Delta_t) Exposure Window & Payload Safing Controller
# + Fast Wireframe / Slow Core Multi-Clock Settlement & Collar Bridge (Turns L22..L28)
# Strictly zero em-dashes and zero en-dashes.
# ==============================================================================

"""
    PayloadSafetyControllerSpec

Configuration contract for automated Health-Physics ALARA payload safing (`PayloadController.jl`),
bounding exposure time `Delta_t` and cumulative radiation dose `D(t)` before biological survival
or metabolic maintenance thresholds are violated.
"""
struct PayloadSafetyControllerSpec
    safing_dose_fraction_threshold::Float64
    min_calibrated_temperature_K::Float64
    max_biocompatible_temperature_K::Float64
    automated_shutter_isolation_enabled::Bool
    require_stoppable_algorithmic_trigger::Bool

    function PayloadSafetyControllerSpec(;
            safing_dose_fraction_threshold::Real = 0.80,
            min_calibrated_temperature_K::Real = 280.0,
            max_biocompatible_temperature_K::Real = 330.0,
            automated_shutter_isolation_enabled::Bool = true,
            require_stoppable_algorithmic_trigger::Bool = true)
        (isfinite(safing_dose_fraction_threshold) && 0.10 <= safing_dose_fraction_threshold < 1.0) ||
            throw(ArgumentError("safing_dose_fraction_threshold must lie in [0.10, 1.0)"))
        (isfinite(min_calibrated_temperature_K) && min_calibrated_temperature_K > 0.0) ||
            throw(ArgumentError("min_calibrated_temperature_K must be > 0"))
        (isfinite(max_biocompatible_temperature_K) && max_biocompatible_temperature_K > min_calibrated_temperature_K) ||
            throw(ArgumentError("max_biocompatible_temperature_K must exceed min_calibrated_temperature_K"))
        automated_shutter_isolation_enabled ||
            throw(ArgumentError("Health-Physics ALARA Violation: automated_shutter_isolation_enabled must be true"))
        new(
            Float64(safing_dose_fraction_threshold),
            Float64(min_calibrated_temperature_K),
            Float64(max_biocompatible_temperature_K),
            automated_shutter_isolation_enabled,
            require_stoppable_algorithmic_trigger
        )
    end
end

"""
    evaluate_alara_payload_controller_safing(
        clock::MissionClock = MissionClock(mission_elapsed_time_s = 3600.0),
        alara_contract::HealthPhysicsALARAContract = HealthPhysicsALARAContract(),
        controller::PayloadSafetyControllerSpec = PayloadSafetyControllerSpec();
        spec::SpectralRadiationFieldSpec = SpectralRadiationFieldSpec(),
        depths_m::AbstractVector{<:Real} = [0.0, 0.002, 0.004, 0.006],
        temperature_K::Real = 298.15,
        spice_units_verified_by_spkpos::Bool = true,
        algorithmic_ladder_bounds::Tuple{<:Real,<:Real} = (900.0, 2500.0),
        algorithmic_exit_trigger_stoppable::Bool = true)

Couples `MissionClock` exposure windows (`Delta_t`) with `HealthPhysicsALARAContract` (`Distance 1/r^2`
and `Shielding exp(-mu*z)` + `Melanin/Neurospora`) to determine automated payload controller state:
- `:EMERGENCY_SCRAM_ISOLATE` (`:Rejected`) if cumulative dose exceeds `max_allowable_cumulative_dose_Gy` or `T > max_biocompatible_temperature_K`.
- `:SAFE_HOLD_SHUTTER_CLOSED` (`:Incomplete`) if cumulative dose reaches `safing_dose_fraction_threshold * max_allowable_cumulative_dose_Gy` (default `0.80 * D_max`).
- `:THERMAL_HIBERNATION_HOLD` (`:Incomplete`) if `temperature_K < min_calibrated_temperature_K` (`280.0 K` NutMEG calibration boundary).
- `:NOMINAL_BIOREACTOR_SCIENCE` (`:Accepted`) when all ALARA, thermal, and SPICE unit gates pass.
Also enforces Turn `L27` (`index = 45`) algorithmic governance: any automated ladder (`900..2500`) must have
an explicit, stoppable exit trigger (`algorithmic_exit_trigger_stoppable == true`).
"""
function evaluate_alara_payload_controller_safing(
        clock::MissionClock = MissionClock(mission_elapsed_time_s = 3600.0),
        alara_contract::HealthPhysicsALARAContract = HealthPhysicsALARAContract(),
        controller::PayloadSafetyControllerSpec = PayloadSafetyControllerSpec();
        spec::SpectralRadiationFieldSpec = SpectralRadiationFieldSpec(),
        depths_m::AbstractVector{<:Real} = [0.0, 0.002, 0.004, 0.006],
        temperature_K::Real = 298.15,
        spice_units_verified_by_spkpos::Bool = true,
        algorithmic_ladder_bounds::Tuple{<:Real,<:Real} = (900.0, 2500.0),
        algorithmic_exit_trigger_stoppable::Bool = true)
    if controller.require_stoppable_algorithmic_trigger && !algorithmic_exit_trigger_stoppable
        throw(ArgumentError(
            "Operational Hippocratism Violation (Turn L27): automated algorithmic ladder ($(algorithmic_ladder_bounds[1])..$(algorithmic_ladder_bounds[2])) cannot run without an explicit stoppable exit-trigger lever"
        ))
    end
    (algorithmic_ladder_bounds[1] > 0.0 && algorithmic_ladder_bounds[2] >= algorithmic_ladder_bounds[1]) ||
        throw(ArgumentError("algorithmic_ladder_bounds must satisfy 0 < base <= max"))

    effective_exposure_s = max(clock.mission_elapsed_time_s, alara_contract.time_exposure_s)
    active_contract = HealthPhysicsALARAContract(
        time_exposure_s = effective_exposure_s,
        max_allowable_cumulative_dose_Gy = alara_contract.max_allowable_cumulative_dose_Gy,
        reference_distance_m = alara_contract.reference_distance_m,
        standoff_distance_m = alara_contract.standoff_distance_m,
        shielding_regolith_thickness_m = alara_contract.shielding_regolith_thickness_m,
        shielding_water_clay_thickness_m = alara_contract.shielding_water_clay_thickness_m,
        melanin_mass_fraction = alara_contract.melanin_mass_fraction,
        neurospora_circadian_gate_active = alara_contract.neurospora_circadian_gate_active,
        epistemic_hippocratism_enforced = alara_contract.epistemic_hippocratism_enforced
    )

    alara_eval = evaluate_health_physics_alara_and_hippocratism(
        spec,
        active_contract;
        depths_m = depths_m,
        temperature_K = temperature_K,
        spice_units_verified_by_spkpos = spice_units_verified_by_spkpos
    )

    peak_dose_Gy = alara_eval.time_control.peak_cumulative_dose_Gy
    d_max_Gy = active_contract.max_allowable_cumulative_dose_Gy
    d_safing_Gy = controller.safing_dose_fraction_threshold * d_max_Gy
    t_K = Float64(temperature_K)

    controller_mode = if (peak_dose_Gy > d_max_Gy) || (t_K > controller.max_biocompatible_temperature_K)
        :EMERGENCY_SCRAM_ISOLATE
    elseif peak_dose_Gy >= d_safing_Gy
        :SAFE_HOLD_SHUTTER_CLOSED
    elseif (t_K < controller.min_calibrated_temperature_K) || !spice_units_verified_by_spkpos
        :THERMAL_HIBERNATION_HOLD
    else
        :NOMINAL_BIOREACTOR_SCIENCE
    end

    shutter_closed = (controller_mode != :NOMINAL_BIOREACTOR_SCIENCE)
    remaining_safe_time_s = max(0.0, (d_safing_Gy - peak_dose_Gy) / max(maximum(alara_eval.dose_rate_Gy_s), 1.0e-30))

    return (
        controller_mode = controller_mode,
        shutter_closed = shutter_closed,
        epistemic_verdict = alara_eval.epistemic_hippocratism_verdict,
        alara_satisfied = alara_eval.alara_satisfied && (controller_mode == :NOMINAL_BIOREACTOR_SCIENCE),
        effective_exposure_time_s = effective_exposure_s,
        peak_cumulative_dose_Gy = peak_dose_Gy,
        safing_threshold_dose_Gy = d_safing_Gy,
        max_allowable_dose_Gy = d_max_Gy,
        remaining_time_to_safing_s = remaining_safe_time_s,
        temperature_K = t_K,
        at_exact_280K_nutmeg_boundary = (t_K == controller.min_calibrated_temperature_K),
        algorithmic_ladder_bounds = (Float64(algorithmic_ladder_bounds[1]), Float64(algorithmic_ladder_bounds[2])),
        algorithmic_exit_trigger_stoppable = algorithmic_exit_trigger_stoppable,
        alara_evaluation = alara_eval
    )
end
