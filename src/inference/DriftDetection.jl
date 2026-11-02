# ==============================================================================
# src/inference/DriftDetection.jl
# Posterior Parameter Mean Drift Delta_t = mu_{theta,t} - mu_{theta,t-1} and
# Covariance Drift Delta_Sigma_t = Sigma_{theta,t} - Sigma_{theta,t-1}
# Stored Strictly Separately from Physical State Change, Priors, and
# Observation-Model Parameters
# ==============================================================================

struct PosteriorDriftRecord
    step_index::Int
    delta_mu_t::Vector{Float64}
    delta_sigma_t::Matrix{Float64}
    mean_l2_drift::Float64
    covariance_frobenius_drift::Float64
    separated_from_physical_state::Bool
    separated_from_prior::Bool
    separated_from_observation_model::Bool
end

"""
    compute_posterior_drift(
        mu_prev::AbstractVector{<:Real},
        sigma_prev::AbstractMatrix{<:Real},
        mu_curr::AbstractVector{<:Real},
        sigma_curr::AbstractMatrix{<:Real};
        step_index::Int = 1,
        overwrite_physical_state_or_prior::Bool = false) -> PosteriorDriftRecord

Computes posterior parameter mean shift `Delta_t = mu_{theta,t} - mu_{theta,t-1}` and
covariance shift `Delta_Sigma_t = Sigma_{theta,t} - Sigma_{theta,t-1}` without mutating
physical solver state, prior specifications, or observation-model calibration.
"""
function compute_posterior_drift(
        mu_prev::AbstractVector{<:Real},
        sigma_prev::AbstractMatrix{<:Real},
        mu_curr::AbstractVector{<:Real},
        sigma_curr::AbstractMatrix{<:Real};
        step_index::Int = 1,
        overwrite_physical_state_or_prior::Bool = false)
    overwrite_physical_state_or_prior && throw(ArgumentError(
        "Drift Storage Separation Violation: Posterior parameter drift (Delta_t, Delta_Sigma_t) must be stored separately from physical state changes, priors, and observation models"
    ))
    step_index >= 1 || throw(ArgumentError("step_index must be >= 1"))
    d = length(mu_prev)
    d >= 1 || throw(ArgumentError("mu_prev cannot be empty"))
    length(mu_curr) == d || throw(ArgumentError("mu_curr length must match mu_prev"))
    (size(sigma_prev) == (d, d) && size(sigma_curr) == (d, d)) ||
        throw(ArgumentError("sigma_prev and sigma_curr must be $(d)x$(d) matrices"))

    delta_mu = Float64.(mu_curr) .- Float64.(mu_prev)
    delta_sigma = Float64.(sigma_curr) .- Float64.(sigma_prev)
    l2_mu = norm(delta_mu)
    frob_sig = norm(delta_sigma)

    return PosteriorDriftRecord(
        step_index,
        delta_mu,
        delta_sigma,
        l2_mu,
        frob_sig,
        true,
        true,
        true
    )
end
