# ==============================================================================
# src/inference/ParetoDesign.jl
# 6-Objective Multi-Fidelity Bioprocess Pareto Design Evaluator (arXiv:2508.10970)
# Objective Vector (all minimized):
#   [-biomass_productivity, metal_toxicity, radiation_damage,
#    energy_consumption, permeability_loss, mass_and_volume]^T
# ==============================================================================

const SIX_PARETO_OBJECTIVE_NAMES = (
    :neg_biomass_productivity,
    :metal_toxicity,
    :radiation_damage,
    :energy_consumption,
    :permeability_loss,
    :mass_and_volume,
)

struct ParetoDesignCandidate
    candidate_id::String
    biomass_productivity::Float64
    metal_toxicity::Float64
    radiation_damage::Float64
    energy_consumption::Float64
    permeability_loss::Float64
    mass_and_volume::Float64
    objective_vector::Vector{Float64}

    function ParetoDesignCandidate(;
            candidate_id::AbstractString,
            biomass_productivity::Real,
            metal_toxicity::Real,
            radiation_damage::Real,
            energy_consumption::Real,
            permeability_loss::Real,
            mass_and_volume::Real)
        isempty(strip(candidate_id)) && throw(ArgumentError("candidate_id cannot be empty"))
        (biomass_productivity >= 0.0 && metal_toxicity >= 0.0 && radiation_damage >= 0.0 &&
         energy_consumption >= 0.0 && permeability_loss >= 0.0 && mass_and_volume > 0.0) ||
            throw(ArgumentError("Invalid physical quantities for ParetoDesignCandidate $(candidate_id)"))

        obj = Float64[
            -Float64(biomass_productivity),
            Float64(metal_toxicity),
            Float64(radiation_damage),
            Float64(energy_consumption),
            Float64(permeability_loss),
            Float64(mass_and_volume)
        ]
        new(
            String(candidate_id),
            Float64(biomass_productivity),
            Float64(metal_toxicity),
            Float64(radiation_damage),
            Float64(energy_consumption),
            Float64(permeability_loss),
            Float64(mass_and_volume),
            obj
        )
    end
end

"""
    pareto_dominates(a::ParetoDesignCandidate, b::ParetoDesignCandidate) -> Bool

Returns `true` iff candidate `a` is no worse than `b` in all 6 minimized objectives
and strictly better in at least one objective.
"""
function pareto_dominates(a::ParetoDesignCandidate, b::ParetoDesignCandidate)
    va = a.objective_vector
    vb = b.objective_vector
    return all(va[i] <= vb[i] for i in 1:6) && any(va[i] < vb[i] for i in 1:6)
end

"""
    extract_pareto_front(candidates::AbstractVector{ParetoDesignCandidate})

Extracts the non-dominated Pareto front across the 6-objective design space.
"""
function extract_pareto_front(candidates::AbstractVector{ParetoDesignCandidate})
    isempty(candidates) && throw(ArgumentError("candidates vector cannot be empty"))
    n = length(candidates)
    is_dominated = falses(n)
    for i in 1:n
        for j in 1:n
            if i != j && pareto_dominates(candidates[j], candidates[i])
                is_dominated[i] = true
                break
            end
        end
    end
    front = [candidates[i] for i in 1:n if !is_dominated[i]]
    return (
        objective_names = SIX_PARETO_OBJECTIVE_NAMES,
        total_candidates = n,
        pareto_front_count = length(front),
        pareto_front = front,
        dominated_candidates = [candidates[i] for i in 1:n if is_dominated[i]]
    )
end
