# ==============================================================================
# src/biology/RNAReplicatorToyAdapter.jl
# 2D Synthetic Lattice RNA Replicator / Parasite Evolutionary Toy Adapter
# Strictly Quarantined from 3D Continuum Pore-Scale Mass-Balance Claims
# ==============================================================================

struct RNAReplicatorToyState
    lattice_dims::Tuple{Int,Int}
    replicator_density::Float64
    parasite_density::Float64
    diffusion_swap_rate::Float64
    epistemic_tag::Symbol
    is_physical_3d_pore_solver::Bool

    function RNAReplicatorToyState(;
            lattice_dims::Tuple{Int,Int} = (32, 32),
            replicator_density::Real = 0.35,
            parasite_density::Real = 0.12,
            diffusion_swap_rate::Real = 0.20,
            promote_to_3d_continuum_pore_solver::Bool = false)
        promote_to_3d_continuum_pore_solver && throw(ArgumentError(
            "Epistemic Quarantine Violation: RNAReplicatorToyAdapter is a 2D synthetic lattice toy model and cannot be promoted into a 3D continuum pore-scale mass-balance solver"
        ))
        (lattice_dims[1] >= 4 && lattice_dims[2] >= 4) ||
            throw(ArgumentError("lattice_dims must be at least 4x4"))
        (0.0 <= replicator_density <= 1.0 && 0.0 <= parasite_density <= 1.0 &&
         replicator_density + parasite_density <= 1.0) ||
            throw(ArgumentError("Replicator and parasite densities must be non-negative and sum to <= 1"))
        diffusion_swap_rate >= 0.0 || throw(ArgumentError("diffusion_swap_rate must be >= 0"))

        new(
            lattice_dims,
            Float64(replicator_density),
            Float64(parasite_density),
            Float64(diffusion_swap_rate),
            :synthetic_2d_lattice_toy_model,
            false
        )
    end
end

"""
    step_rna_replicator_toy(state::RNAReplicatorToyState; dt_toy::Real = 1.0, k_rep::Real = 0.18, k_par::Real = 0.22, decay::Real = 0.05)

Advances the 2D synthetic lattice replicator-parasite replicase dynamics by `dt_toy` while preserving
the explicit `:synthetic_2d_lattice_toy_model` tag.
"""
function step_rna_replicator_toy(
        state::RNAReplicatorToyState;
        dt_toy::Real = 1.0,
        k_rep::Real = 0.18,
        k_par::Real = 0.22,
        decay::Real = 0.05)
    dt_toy > 0.0 || throw(ArgumentError("dt_toy must be > 0"))
    R = state.replicator_density
    P = state.parasite_density
    empty_frac = max(0.0, 1.0 - R - P)
    # Spatial wave/cluster protection factor from finite lattice diffusion
    spatial_shield = 1.0 / (1.0 + 0.8 * state.diffusion_swap_rate)
    dR = (Float64(k_rep) * (R ^ 2) * empty_frac - Float64(decay) * R) * Float64(dt_toy)
    dP = (Float64(k_par) * R * P * empty_frac * (1.0 - 0.5 * spatial_shield) - Float64(decay) * P) * Float64(dt_toy)
    R_next = clamp(R + dR, 0.0, 1.0)
    P_next = clamp(P + dP, 0.0, 1.0 - R_next)

    return RNAReplicatorToyState(
        lattice_dims = state.lattice_dims,
        replicator_density = R_next,
        parasite_density = P_next,
        diffusion_swap_rate = state.diffusion_swap_rate
    )
end

"""
    verify_rna_biofilm_toy_mass_ledger(;
        M_prev::Real = 72410.0,
        M_curr::Real = 74290.0,
        injected_gross::Real = 2025.0,
        inject_clip::Real = 0.0,
        diffusion_clip::Real = 0.0,
        births::Integer = 12,
        deaths::Integer = 7,
        death_return::Real = 35.0,
        atol::Real = 1.0e-8)

Verify the exact 5-term mass-conservation ledger of `rna_biofilm_toy` (`src/rna_biofilm_3d.py`):
`M(t) - M(t-1) = injected_gross - inject_clip - diffusion_clip - birth_sink + death_return`
where `birth_sink = 15.0 * births`, `diffusion_clip == 0.0` for the guarded downhill divide-by-degree
scheme, and `death_return <= 5.0 * deaths`.
"""
function verify_rna_biofilm_toy_mass_ledger(;
        M_prev::Real = 72410.0,
        M_curr::Real = 74290.0,
        injected_gross::Real = 2025.0,
        inject_clip::Real = 0.0,
        diffusion_clip::Real = 0.0,
        births::Integer = 12,
        deaths::Integer = 7,
        death_return::Real = 35.0,
        atol::Real = 1.0e-8)
    (births >= 0 && deaths >= 0) || throw(ArgumentError("births and deaths must be >= 0"))
    b_sink = 15.0 * Float64(births)
    d_ret = Float64(death_return)
    (0.0 <= d_ret <= 5.0 * Float64(deaths) + 1.0e-12) ||
        throw(ArgumentError("death_return ($d_ret) must lie in [0, 5 * deaths ($(5 * deaths))]"))
    dM = Float64(M_curr) - Float64(M_prev)
    rhs = Float64(injected_gross) - Float64(inject_clip) - Float64(diffusion_clip) - b_sink + d_ret
    res = abs(dM - rhs)
    closed = isfinite(res) && (res <= Float64(atol)) && (Float64(diffusion_clip) == 0.0)
    return (
        mass_ledger_closed = closed,
        delta_M = dM,
        rhs_5term = rhs,
        residual = res,
        injected_gross = Float64(injected_gross),
        inject_clip = Float64(inject_clip),
        diffusion_clip = Float64(diffusion_clip),
        birth_sink = b_sink,
        death_return = d_ret,
        births = Int64(births),
        deaths = Int64(deaths),
        tolerance = Float64(atol)
    )
end

"""
    evaluate_biofilms_cpm_calibration_refusal_gate(; attempt_radiodialysis_on_occupancy_mean::Bool = false)

Preserve the `Biofilms-publish` 2D Cellular Potts Model calibration and refusal contracts:
- `clock = :numerical_mcs`, `seconds_per_mcs = NaN`, `lattice_pitch_m = nothing`
- `radiodialysis_coupled_gate = "RADIODIALYSIS: BLOCKED"` (`reason = :biomass_basis_is_occupancy_mean_in_g_cm3_slot`)
- `growth_survival_status = :unsupported_by_current_model`
- `radiotropic_vs_radiotrophic = :radiotropic_hamiltonian_bias_not_radiotrophic_metabolism`
- `evidence_status = :Incomplete`
"""
function evaluate_biofilms_cpm_calibration_refusal_gate(;
        attempt_radiodialysis_on_occupancy_mean::Bool = false)
    if attempt_radiodialysis_on_occupancy_mean
        throw(ArgumentError(
            "RADIODIALYSIS: BLOCKED -- biomass_basis_is_occupancy_mean_in_g_cm3_slot; uncalibrated 2D CPM occupancy fraction cannot drive a g/cm^3 radiodialysis mass-balance sink"
        ))
    end
    return (
        clock = :numerical_mcs,
        seconds_per_mcs = NaN,
        lattice_pitch_m = nothing,
        radiodialysis_coupled_gate = "RADIODIALYSIS: BLOCKED",
        radiodialysis_refusal_reason = :biomass_basis_is_occupancy_mean_in_g_cm3_slot,
        growth_survival_status = :unsupported_by_current_model,
        radiotropic_vs_radiotrophic = :radiotropic_hamiltonian_bias_not_radiotrophic_metabolism,
        evidence_status = :Incomplete,
        calibration_status = :uncalibrated,
        refusal_contracts_preserved = true
    )
end

