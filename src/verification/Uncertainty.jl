# ==============================================================================
# src/verification/Uncertainty.jl
# Epistemic vs Aleatoric Uncertainty Decomposition, Equifinality Audit, and
# Multivariate Out-of-Distribution (OOD) Shift Gating
# ==============================================================================

"""
    decompose_epistemic_aleatoric_uncertainty(
        ensemble_predictions::AbstractMatrix{<:Real},
        observation_noise_variance::Real)

Decomposes total predictive variance across an ensemble of parameter draws (`N_obs x N_draws`)
into epistemic variance (between-draw model spread) and aleatoric variance (measurement noise).
"""
function decompose_epistemic_aleatoric_uncertainty(
        ensemble_predictions::AbstractMatrix{<:Real},
        observation_noise_variance::Real)
    (size(ensemble_predictions, 1) >= 1 && size(ensemble_predictions, 2) >= 2) ||
        throw(ArgumentError("ensemble_predictions must have size (N_obs >= 1, N_draws >= 2)"))
    observation_noise_variance >= 0.0 ||
        throw(ArgumentError("observation_noise_variance must be >= 0"))

    epistemic_var = vec(var(Float64.(ensemble_predictions); dims = 2))
    aleatoric_var = fill(Float64(observation_noise_variance), length(epistemic_var))
    total_var = epistemic_var .+ aleatoric_var

    return (
        epistemic_variance = epistemic_var,
        aleatoric_variance = aleatoric_var,
        total_predictive_variance = total_var,
        epistemic_fraction = [total_var[i] > 0.0 ? epistemic_var[i] / total_var[i] : 0.0 for i in eachindex(total_var)]
    )
end

"""
    audit_equifinality_underdetermination(
        theta_1::AbstractVector{<:Real},
        theta_2::AbstractVector{<:Real},
        forward_fn::Function;
        obs_tol::Real = 1.0e-9)

Tests whether two distinct parameter vectors `theta_1 != theta_2` map to indistinguishable
terminal observables `forward_fn(theta_1) == forward_fn(theta_2)` (equifinality), forbidding
unique backward parameter inversion without additional orthogonal measurements.
"""
function audit_equifinality_underdetermination(
        theta_1::AbstractVector{<:Real},
        theta_2::AbstractVector{<:Real},
        forward_fn::Function;
        obs_tol::Real = 1.0e-9)
    length(theta_1) == length(theta_2) || throw(ArgumentError("theta_1 and theta_2 must have equal length"))
    obs_tol > 0.0 || throw(ArgumentError("obs_tol must be > 0"))

    param_dist = norm(Float64.(theta_1) .- Float64.(theta_2))
    y1 = Float64.(forward_fn(Float64.(theta_1)))
    y2 = Float64.(forward_fn(Float64.(theta_2)))
    obs_diff = norm(y1 .- y2)
    is_equifinal = (param_dist > 1.0e-6) && (obs_diff <= Float64(obs_tol))

    return (
        parameters_distinct = (param_dist > 1.0e-6),
        parameter_l2_distance = param_dist,
        observable_l2_difference = obs_diff,
        is_equifinal = is_equifinal,
        unique_backward_inversion_forbidden = is_equifinal
    )
end
