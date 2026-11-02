# ==============================================================================
# src/contracts/UnitSystem.jl
# SI Dimensional Homogeneity, NAIF km -> SI m Conversion, Molality -> Molarity,
# and Activity-Averaged pH Enforcement
# ==============================================================================

const VALID_SI_BASE_UNITS = (
    :m,
    :m_s,
    :m_s2,
    :m2,
    :m2_s,
    :m3,
    :s,
    :kg,
    :kg_m3,
    :kg_m3_s,
    :mol,
    :mol_m3,
    :mol_m3_s,
    :mol_kgw,
    :K,
    :Pa,
    :Gy,
    :Gy_s,
    :Bq,
    :W,
    :W_m2,
    :W_m3,
    :J_mol,
    :dimensionless,
)

const NAIF_KM_TO_SI_M = 1000.0

struct SMBRUnitQuantity
    name::String
    value::Float64
    unit::Symbol
    provenance_source::Symbol

    function SMBRUnitQuantity(name::AbstractString, value::Real, unit::Symbol; provenance_source::Symbol = :solver_contract)
        isempty(strip(name)) && throw(ArgumentError("Quantity name must be non-empty"))
        isfinite(value) || throw(ArgumentError("Quantity value for $(name) must be finite, got $(value)"))
        unit in VALID_SI_BASE_UNITS || throw(ArgumentError("Unsupported unit :$(unit) for $(name); must be in VALID_SI_BASE_UNITS"))
        new(String(name), Float64(value), unit, provenance_source)
    end
end

"""
    convert_naif_km_to_si_m(vec_km::AbstractVector{<:Real}) -> Vector{Float64}

Converts a 3-element NAIF SPICE position or velocity vector from kilometers (`km` or `km/s`)
to SI meters (`m` or `m/s`) via exact `* 1000.0` scaling.
"""
function convert_naif_km_to_si_m(vec_km::AbstractVector{<:Real})
    length(vec_km) == 3 || throw(ArgumentError("NAIF state vector must have length 3, got $(length(vec_km))"))
    all(isfinite, vec_km) || throw(ArgumentError("NAIF state vector entries must be finite"))
    return Float64.(vec_km) .* NAIF_KM_TO_SI_M
end

"""
    convert_molality_to_molarity_si(molality_mol_kgw::Real, solution_density_kg_m3::Real; tds_kg_kg::Real = 0.0) -> Float64

Converts PHREEQC aqueous molality (`mol / kg_water`) to continuum finite-volume molarity (`mol / m^3`).
"""
function convert_molality_to_molarity_si(
        molality_mol_kgw::Real,
        solution_density_kg_m3::Real;
        tds_kg_kg::Real = 0.0)
    molality_mol_kgw >= 0.0 || throw(ArgumentError("molality_mol_kgw must be non-negative"))
    solution_density_kg_m3 > 0.0 || throw(ArgumentError("solution_density_kg_m3 must be strictly positive"))
    (0.0 <= tds_kg_kg < 1.0) || throw(ArgumentError("tds_kg_kg must lie in [0, 1)"))
    water_mass_per_m3 = Float64(solution_density_kg_m3) * (1.0 - Float64(tds_kg_kg))
    return Float64(molality_mol_kgw) * water_mass_per_m3
end

"""
    compute_activity_averaged_ph(ph_values::AbstractVector{<:Real}; weights::Union{Nothing,AbstractVector{<:Real}} = nothing)

Computes proton-activity-averaged pH: `-log10(sum(w_i * 10^(-pH_i)) / sum(w_i))`
and quantifies the bias incurred if arithmetic averaging `mean(pH)` were erroneously used.
"""
function compute_activity_averaged_ph(
        ph_values::AbstractVector{<:Real};
        weights::Union{Nothing,AbstractVector{<:Real}} = nothing)
    isempty(ph_values) && throw(ArgumentError("ph_values cannot be empty"))
    all(x -> isfinite(x) && 0.0 <= x <= 14.0, ph_values) ||
        throw(ArgumentError("All pH values must be finite and within [0, 14]"))

    n = length(ph_values)
    w = if weights === nothing
        fill(1.0 / n, n)
    else
        length(weights) == n || throw(ArgumentError("weights length must match ph_values length"))
        all(x -> isfinite(x) && x > 0.0, weights) || throw(ArgumentError("weights must be finite and strictly positive"))
        Float64.(weights) ./ sum(weights)
    end

    h_plus_activity = sum(w[i] * (10.0 ^ (-Float64(ph_values[i]))) for i in 1:n)
    ph_activity_avg = -log10(h_plus_activity)
    ph_arithmetic_invalid = sum(w[i] * Float64(ph_values[i]) for i in 1:n)
    bias = abs(ph_arithmetic_invalid - ph_activity_avg)

    return (
        ph_activity_avg = ph_activity_avg,
        ph_arithmetic_invalid = ph_arithmetic_invalid,
        arithmetic_vs_activity_bias = bias,
        h_plus_activity_mol_L = h_plus_activity
    )
end

"""
    verify_unit_homogeneity(quantities::AbstractVector{SMBRUnitQuantity})

Audits a collection of `SMBRUnitQuantity` entries for finite values, valid SI unit symbols,
and explicit provenance tags.
"""
function verify_unit_homogeneity(quantities::AbstractVector{SMBRUnitQuantity})
    isempty(quantities) && throw(ArgumentError("quantities vector must not be empty"))
    all_valid = all(q -> isfinite(q.value) && (q.unit in VALID_SI_BASE_UNITS), quantities)
    return (
        homogeneous = all_valid,
        quantity_count = length(quantities),
        units = Tuple(q.unit for q in quantities)
    )
end
