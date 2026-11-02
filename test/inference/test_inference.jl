@testset "Inference: BayesianCalibration, HastingsMetropolis, DriftDetection, ParetoDesign" begin
    prior = GaussianPriorSpec()
    lp = evaluate_log_posterior([0.12, 25.0, 1.5e-6], [0.12, 25.0, 1.5e-6], identity, prior)
    @test isfinite(lp.log_posterior_unnormalized)

    alpha_1 = compute_metropolis_hastings_alpha(-10.0, -8.0)
    alpha_2 = compute_metropolis_hastings_alpha(-8.0, -10.0)
    @test alpha_1 == 1.0
    @test isapprox(alpha_2, exp(-2.0); atol = 1.0e-12)

    mcmc = run_hastings_metropolis_calibration(n_steps = 200, seed = 101)
    @test mcmc.sampling_acceptance_is_not_engineering_verdict
    @test 0.0 < mcmc.acceptance_rate <= 1.0
    @test mcmc.min_ess > 0.0
    @test isfinite(mcmc.ppc_rmse)
    @test_throws ArgumentError run_hastings_metropolis_calibration(
        conflate_mcmc_acceptance_with_engineering_build_verdict = true
    )

    drift = compute_posterior_drift(
        [0.12, 25.0, 1.5e-6],
        Matrix{Float64}(I, 3, 3) .* 0.01,
        mcmc.posterior_mean,
        mcmc.posterior_cov
    )
    @test drift.separated_from_physical_state
    @test drift.separated_from_prior
    @test drift.separated_from_observation_model
    @test drift.mean_l2_drift >= 0.0
    @test_throws ArgumentError compute_posterior_drift(
        [0.12, 25.0, 1.5e-6],
        Matrix{Float64}(I, 3, 3),
        [0.12, 25.0, 1.5e-6],
        Matrix{Float64}(I, 3, 3);
        overwrite_physical_state_or_prior = true
    )

    c1 = ParetoDesignCandidate(
        candidate_id = "c1_pareto_optimal",
        biomass_productivity = 0.90,
        metal_toxicity = 0.04,
        radiation_damage = 0.03,
        energy_consumption = 85.0,
        permeability_loss = 0.10,
        mass_and_volume = 3.8
    )
    c2 = ParetoDesignCandidate(
        candidate_id = "c2_pareto_optimal_low_energy",
        biomass_productivity = 0.70,
        metal_toxicity = 0.04,
        radiation_damage = 0.03,
        energy_consumption = 55.0,
        permeability_loss = 0.10,
        mass_and_volume = 3.5
    )
    c3 = ParetoDesignCandidate(
        candidate_id = "c3_dominated",
        biomass_productivity = 0.65,
        metal_toxicity = 0.10,
        radiation_damage = 0.08,
        energy_consumption = 95.0,
        permeability_loss = 0.25,
        mass_and_volume = 5.0
    )
    @test pareto_dominates(c1, c3)
    @test !pareto_dominates(c1, c2)
    pf = extract_pareto_front([c1, c2, c3])
    @test pf.pareto_front_count == 2
    @test length(pf.dominated_candidates) == 1
end
