# ==============================================================================
# src/inference/BayesianCalibration.jl
# Bayesian Posterior Calibration p(theta | y) propto p(y | theta) * p(theta)
# and Multi-Fidelity Batch Bayesian Optimization (arXiv:2508.10970)
# ==============================================================================

struct GaussianPriorSpec
    param_names::Vector{Symbol}
    prior_mean::Vector{Float64}
    prior_std::Vector{Float64}

    function GaussianPriorSpec(;
            param_names::AbstractVector{Symbol} = [:mu_max, :K_bind, :chi_R],
            prior_mean::AbstractVector{<:Real} = [0.12, 25.0, 1.5e-6],
            prior_std::AbstractVector{<:Real} = [0.04, 8.0, 0.5e-6])
        n = length(param_names)
        n >= 1 || throw(ArgumentError("param_names cannot be empty"))
        (length(prior_mean) == n && length(prior_std) == n) ||
            throw(ArgumentError("prior_mean and prior_std must match param_names length"))
        all(s -> isfinite(s) && s > 0.0, prior_std) ||
            throw(ArgumentError("prior_std entries must be finite and > 0"))
        new(Vector{Symbol}(param_names), Float64.(prior_mean), Float64.(prior_std))
    end
end

"""
    evaluate_log_posterior(theta::AbstractVector{<:Real}, y_obs::AbstractVector{<:Real},
        forward_pred_fn::Function, prior::GaussianPriorSpec; obs_noise_std::Real = 0.05)

Computes `log p(theta | y) = log p(y | theta) + log p(theta)` (up to normalizing constant `log p(y)`).
"""
function evaluate_log_posterior(
        theta::AbstractVector{<:Real},
        y_obs::AbstractVector{<:Real},
        forward_pred_fn::Function,
        prior::GaussianPriorSpec;
        obs_noise_std::Real = 0.05)
    length(theta) == length(prior.param_names) ||
        throw(ArgumentError("theta length must match prior.param_names"))
    isempty(y_obs) && throw(ArgumentError("y_obs cannot be empty"))
    obs_noise_std > 0.0 || throw(ArgumentError("obs_noise_std must be > 0"))

    log_prior = 0.0
    for i in eachindex(theta)
        z = (Float64(theta[i]) - prior.prior_mean[i]) / prior.prior_std[i]
        log_prior += -0.5 * (z ^ 2) - log(sqrt(2.0 * pi) * prior.prior_std[i])
    end

    y_pred = forward_pred_fn(Float64.(theta))
    length(y_pred) == length(y_obs) ||
        throw(ArgumentError("forward_pred_fn output length must match y_obs length"))

    sig = Float64(obs_noise_std)
    log_lik = 0.0
    for j in eachindex(y_obs)
        r = (Float64(y_obs[j]) - Float64(y_pred[j])) / sig
        log_lik += -0.5 * (r ^ 2) - log(sqrt(2.0 * pi) * sig)
    end

    return (
        log_prior = log_prior,
        log_likelihood = log_lik,
        log_posterior_unnormalized = log_prior + log_lik,
        y_predicted = Float64.(y_pred)
    )
end
