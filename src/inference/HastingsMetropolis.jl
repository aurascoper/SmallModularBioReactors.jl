# ==============================================================================
# src/inference/HastingsMetropolis.jl
# Metropolis-Hastings MCMC Sampler with Acceptance Ratio alpha(theta, theta'),
# Effective Sample Size (ESS), Acceptance Rate, Lag-1 Autocorrelation, and
# Posterior Predictive Check (PPC)
# Strictly Separated from Engineering Build / Deployment Decisions
# ==============================================================================

"""
    compute_metropolis_hastings_alpha(log_post_current::Real, log_post_proposed::Real;
        log_q_current_given_proposed::Real = 0.0,
        log_q_proposed_given_current::Real = 0.0) -> Float64

Computes the Hastings-Metropolis acceptance probability:
`alpha(theta, theta') = min(1, exp((log_post_proposed + log_q(theta | theta')) - (log_post_current + log_q(theta' | theta))))`
"""
function compute_metropolis_hastings_alpha(
        log_post_current::Real,
        log_post_proposed::Real;
        log_q_current_given_proposed::Real = 0.0,
        log_q_proposed_given_current::Real = 0.0)
    (isfinite(log_post_current) && isfinite(log_post_proposed)) ||
        throw(ArgumentError("Log-posterior values must be finite"))
    log_ratio = (Float64(log_post_proposed) + Float64(log_q_current_given_proposed)) -
                (Float64(log_post_current) + Float64(log_q_proposed_given_current))
    return log_ratio >= 0.0 ? 1.0 : exp(log_ratio)
end

"""
    run_hastings_metropolis_calibration(;
        y_obs::AbstractVector{<:Real} = [0.24, 24.5, 1.45e-6],
        prior::GaussianPriorSpec = GaussianPriorSpec(),
        forward_pred_fn::Function = identity,
        obs_noise_std::Real = 0.05,
        proposal_scales::AbstractVector{<:Real} = [0.015, 2.5, 0.2e-6],
        n_steps::Int = 400,
        seed::Int = 42,
        conflate_mcmc_acceptance_with_engineering_build_verdict::Bool = false)

Runs a reproducible Random-Walk Metropolis-Hastings chain and returns posterior samples,
mean, covariance, acceptance rate, lag-1 autocorrelation, Effective Sample Size (`ESS`),
and Posterior Predictive Check (`PPC` RMSE).
Throws `ArgumentError` if `conflate_mcmc_acceptance_with_engineering_build_verdict == true`.
"""
function run_hastings_metropolis_calibration(;
        y_obs::AbstractVector{<:Real} = [0.13, 24.5, 1.45e-6],
        prior::GaussianPriorSpec = GaussianPriorSpec(),
        forward_pred_fn::Function = identity,
        obs_noise_std::Real = 0.05,
        proposal_scales::AbstractVector{<:Real} = [0.015, 2.5, 0.2e-6],
        n_steps::Int = 400,
        seed::Int = 42,
        conflate_mcmc_acceptance_with_engineering_build_verdict::Bool = false)
    conflate_mcmc_acceptance_with_engineering_build_verdict && throw(ArgumentError(
        "Architectural Separation Violation: Hastings-Metropolis proposal acceptance alpha(theta, theta') is a Markov-chain transition kernel, NOT an engineering build or deployment verdict"
    ))
    n_steps >= 20 || throw(ArgumentError("n_steps must be >= 20"))
    d = length(prior.param_names)
    length(proposal_scales) == d || throw(ArgumentError("proposal_scales length must match parameter dimension"))
    all(s -> s > 0.0, proposal_scales) || throw(ArgumentError("proposal_scales must be > 0"))

    rng = MersenneTwister(seed)
    theta_curr = copy(prior.prior_mean)
    eval_curr = evaluate_log_posterior(theta_curr, y_obs, forward_pred_fn, prior; obs_noise_std = obs_noise_std)

    chain = zeros(Float64, d, n_steps)
    accepted_count = 0
    for step in 1:n_steps
        theta_prop = [theta_curr[j] + Float64(proposal_scales[j]) * randn(rng) for j in 1:d]
        eval_prop = evaluate_log_posterior(theta_prop, y_obs, forward_pred_fn, prior; obs_noise_std = obs_noise_std)
        alpha = compute_metropolis_hastings_alpha(
            eval_curr.log_posterior_unnormalized,
            eval_prop.log_posterior_unnormalized
        )
        if rand(rng) < alpha
            theta_curr = theta_prop
            eval_curr = eval_prop
            accepted_count += 1
        end
        chain[:, step] .= theta_curr
    end

    acceptance_rate = accepted_count / n_steps
    burn_in = max(1, div(n_steps, 4))
    post_samples = chain[:, (burn_in + 1):end]
    n_post = size(post_samples, 2)

    posterior_mean = vec(mean(post_samples; dims = 2))
    centered = post_samples .- posterior_mean
    posterior_cov = (centered * centered') ./ max(1, n_post - 1)

    # Lag-1 autocorrelation & ESS per dimension
    lag1_autocorr = zeros(Float64, d)
    ess_per_dim = zeros(Float64, d)
    for j in 1:d
        x = view(post_samples, j, :)
        v = var(x)
        if v <= 1.0e-24
            lag1_autocorr[j] = 0.0
            ess_per_dim[j] = Float64(n_post)
        else
            m = posterior_mean[j]
            num = sum((x[t] - m) * (x[t - 1] - m) for t in 2:n_post) / (n_post - 1)
            rho1 = clamp(num / v, -0.95, 0.98)
            lag1_autocorr[j] = rho1
            ess_per_dim[j] = Float64(n_post) * (1.0 - rho1) / (1.0 + rho1)
        end
    end

    y_ppc = forward_pred_fn(posterior_mean)
    ppc_rmse = sqrt(mean((Float64.(y_obs) .- Float64.(y_ppc)) .^ 2))

    return (
        posterior_mean = posterior_mean,
        posterior_cov = posterior_cov,
        acceptance_rate = acceptance_rate,
        lag1_autocorrelation = lag1_autocorr,
        ess_per_dim = ess_per_dim,
        min_ess = minimum(ess_per_dim),
        ppc_predicted = y_ppc,
        ppc_rmse = ppc_rmse,
        sampling_acceptance_is_not_engineering_verdict = true
    )
end
