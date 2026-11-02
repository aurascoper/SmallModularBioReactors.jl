# ==============================================================================
# src/transport/PHREEQCAdapter.jl
# PHREEQC / PhreeqcRM Cell-Scale Geochemical Speciation, Charge Balance,
# Mineral Saturation Index (SI), and Activity-Averaged pH Adapter
# ==============================================================================

struct PHREEQCCellSpeciation
    cell_index::Int
    pH::Float64
    pe::Float64
    temperature_K::Float64
    cation_equivalents_mol_m3::Float64
    anion_equivalents_mol_m3::Float64
    calcite_SI::Float64
    cerium_oxalate_SI::Float64
    surface_OCa_plus_mol_m3::Float64
    surface_O_minus_mol_m3::Float64

    function PHREEQCCellSpeciation(;
            cell_index::Int = 1,
            pH::Real = 7.15,
            pe::Real = 3.25,
            temperature_K::Real = 298.15,
            cation_equivalents_mol_m3::Real = 4.20,
            anion_equivalents_mol_m3::Real = 4.20,
            calcite_SI::Real = 0.05,
            cerium_oxalate_SI::Real = 0.18,
            surface_OCa_plus_mol_m3::Real = 0.32,
            surface_O_minus_mol_m3::Real = 0.14)
        cell_index >= 1 || throw(ArgumentError("cell_index must be >= 1"))
        (0.0 <= pH <= 14.0) || throw(ArgumentError("pH must lie in [0, 14]"))
        temperature_K > 0.0 || throw(ArgumentError("temperature_K must be > 0"))
        (cation_equivalents_mol_m3 > 0.0 && anion_equivalents_mol_m3 > 0.0) ||
            throw(ArgumentError("Ion equivalents must be > 0"))
        (surface_OCa_plus_mol_m3 >= 0.0 && surface_O_minus_mol_m3 >= 0.0) ||
            throw(ArgumentError("Surface complexation concentrations must be >= 0"))

        new(
            cell_index,
            Float64(pH),
            Float64(pe),
            Float64(temperature_K),
            Float64(cation_equivalents_mol_m3),
            Float64(anion_equivalents_mol_m3),
            Float64(calcite_SI),
            Float64(cerium_oxalate_SI),
            Float64(surface_OCa_plus_mol_m3),
            Float64(surface_O_minus_mol_m3)
        )
    end
end

"""
    evaluate_phreeqc_cell_array(cells::AbstractVector{PHREEQCCellSpeciation}; max_charge_imbalance_pct::Real = 1.0)

Audits an array of PHREEQC cell speciation states for electrical charge balance,
activity-averaged pH (`-log10(mean(10^(-pH)))`), and surface complexation site conservation.
"""
function evaluate_phreeqc_cell_array(
        cells::AbstractVector{PHREEQCCellSpeciation};
        max_charge_imbalance_pct::Real = 1.0)
    isempty(cells) && throw(ArgumentError("cells cannot be empty"))
    max_charge_imbalance_pct > 0.0 || throw(ArgumentError("max_charge_imbalance_pct must be > 0"))

    imbalances_pct = [
        100.0 * abs(c.cation_equivalents_mol_m3 - c.anion_equivalents_mol_m3) /
        (c.cation_equivalents_mol_m3 + c.anion_equivalents_mol_m3)
        for c in cells
    ]
    max_imbalance = maximum(imbalances_pct)
    ph_audit = compute_activity_averaged_ph([c.pH for c in cells])
    total_surface_sites = sum(c.surface_OCa_plus_mol_m3 + c.surface_O_minus_mol_m3 for c in cells)

    return (
        cell_count = length(cells),
        charge_balanced = (max_imbalance <= Float64(max_charge_imbalance_pct)),
        max_charge_imbalance_pct = max_imbalance,
        ph_activity_avg = ph_audit.ph_activity_avg,
        ph_arithmetic_invalid = ph_audit.ph_arithmetic_invalid,
        arithmetic_vs_activity_bias = ph_audit.arithmetic_vs_activity_bias,
        total_surface_sites_mol_m3 = total_surface_sites
    )
end
