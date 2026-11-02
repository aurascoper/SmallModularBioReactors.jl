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
