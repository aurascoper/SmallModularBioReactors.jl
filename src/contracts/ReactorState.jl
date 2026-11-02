# ==============================================================================
# src/contracts/ReactorState.jl
# Coupled State X = [u, p, T, phi, k, C_aq, M_min, N_iso, B, E, P],
# Global Conservation dM_i^total/dt = J_i,in - J_i,out + R_i^nuclear,
# 3-Fraction ICP-MS + XRD/XRF + Flow Cytometry Accounting, and System A/B/C Clay
# ==============================================================================

const VALID_CLAY_SYSTEMS = (
    :system_a_sudoite_illite,
    :system_b_bentonite_montmorillonite,
    :system_c_illite_smectite_eps_biofilm,
)

const VALID_ANALYTICAL_INSTRUMENTS = (
    :ICP_MS,
    :XRD_XRF,
    :Flow_Cytometry,
)

struct ThreeFractionAnalyticalLedger
    element_symbol::Symbol
    aqueous_dissolved_mol::Float64
    mineral_adsorbed_mol::Float64
    cellular_polymer_mol::Float64
    icp_ms_instrument::Symbol
    xrd_xrf_instrument::Symbol
    flow_cytometry_instrument::Symbol

    function ThreeFractionAnalyticalLedger(;
            element_symbol::Symbol = :Pb,
            aqueous_dissolved_mol::Real = 1.0e-6,
            mineral_adsorbed_mol::Real = 2.5e-6,
            cellular_polymer_mol::Real = 1.5e-6)
        (aqueous_dissolved_mol >= 0.0 && mineral_adsorbed_mol >= 0.0 && cellular_polymer_mol >= 0.0) ||
            throw(ArgumentError("Analytical fraction molar inventories must be non-negative"))
        new(
            element_symbol,
            Float64(aqueous_dissolved_mol),
            Float64(mineral_adsorbed_mol),
            Float64(cellular_polymer_mol),
            :ICP_MS,
            :XRD_XRF,
            :Flow_Cytometry
        )
    end
end

total_element_moles(ledger::ThreeFractionAnalyticalLedger) =
    ledger.aqueous_dissolved_mol + ledger.mineral_adsorbed_mol + ledger.cellular_polymer_mol

"""
    ReactorState

Coupled state vector `X = [u, p, T, phi, k, C_aq, M_min, N_iso, B, E, P]` across `N` finite-volume cells:
- `u_m_s`: Pore fluid velocity vector per cell (`3 x N` matrix, `m/s`)
- `p_Pa`: Fluid pressure (`N` vector, `Pa`)
- `T_K`: Temperature (`N` vector, `K`)
- `phi`: Porosity (`N` vector, dimensionless in `(0, 1)`)
- `k_m2`: Permeability (`N` vector, `m^2`)
- `C_aq_mol_m3`: Dissolved aqueous species concentrations (`Dict{Symbol,Vector{Float64}}`, `mol/m^3` pore water)
- `M_min_mol_m3`: Cell-scale mineral + surface-adsorbed inventories (`Dict{Symbol,Vector{Float64}}`, `mol/m^3` bulk cell)
- `N_iso_mol_m3`: Isotope / radiolytic species inventories (`Dict{Symbol,Vector{Float64}}`, `mol/m^3` bulk cell)
- `B_kg_m3`: Active microbial biomass density (`N` vector, `kg/m^3`)
- `E_kg_m3`: Extracellular polymeric substance (EPS) density (`N` vector, `kg/m^3`)
- `P_mol_m3`: Metal-binding protein / peptide ligand concentration (`N` vector, `mol/m^3`)
- `clay_system`: One of `:system_a_sudoite_illite`, `:system_b_bentonite_montmorillonite`, `:system_c_illite_smectite_eps_biofilm`
"""
struct ReactorState
    cell_count::Int
    u_m_s::Matrix{Float64}
    p_Pa::Vector{Float64}
    T_K::Vector{Float64}
    phi::Vector{Float64}
    k_m2::Vector{Float64}
    C_aq_mol_m3::Dict{Symbol,Vector{Float64}}
    M_min_mol_m3::Dict{Symbol,Vector{Float64}}
    N_iso_mol_m3::Dict{Symbol,Vector{Float64}}
    B_kg_m3::Vector{Float64}
    E_kg_m3::Vector{Float64}
    P_mol_m3::Vector{Float64}
    clay_system::Symbol

    function ReactorState(;
            cell_count::Int = 4,
            u_m_s::AbstractMatrix{<:Real} = repeat([1.0e-4; 0.0; 0.0], 1, cell_count),
            p_Pa::AbstractVector{<:Real} = fill(101325.0, cell_count),
            T_K::AbstractVector{<:Real} = fill(298.15, cell_count),
            phi::AbstractVector{<:Real} = fill(0.38, cell_count),
            k_m2::AbstractVector{<:Real} = fill(2.5e-13, cell_count),
            C_aq_mol_m3::Dict{Symbol,Vector{Float64}} = Dict(
                :Pb2_plus => fill(0.05, cell_count),
                :Ca2_plus => fill(1.20, cell_count),
                :SO4_2_minus => fill(0.80, cell_count),
                :H2 => fill(0.15, cell_count)
            ),
            M_min_mol_m3::Dict{Symbol,Vector{Float64}} = Dict(
                :Pb2_plus => fill(0.10, cell_count),
                :Ca2_plus => fill(4.50, cell_count)
            ),
            N_iso_mol_m3::Dict{Symbol,Vector{Float64}} = Dict(
                :H2_radiolytic => fill(0.01, cell_count),
                :H2O2_radiolytic => fill(0.002, cell_count)
            ),
            B_kg_m3::AbstractVector{<:Real} = fill(0.45, cell_count),
            E_kg_m3::AbstractVector{<:Real} = fill(0.22, cell_count),
            P_mol_m3::AbstractVector{<:Real} = fill(0.04, cell_count),
            clay_system::Symbol = :system_c_illite_smectite_eps_biofilm)
        cell_count >= 1 || throw(ArgumentError("cell_count must be >= 1"))
        size(u_m_s) == (3, cell_count) || throw(ArgumentError("u_m_s must be 3 x cell_count"))
        length(p_Pa) == cell_count || throw(ArgumentError("p_Pa length mismatch"))
        length(T_K) == cell_count || throw(ArgumentError("T_K length mismatch"))
        length(phi) == cell_count || throw(ArgumentError("phi length mismatch"))
        length(k_m2) == cell_count || throw(ArgumentError("k_m2 length mismatch"))
        length(B_kg_m3) == cell_count || throw(ArgumentError("B_kg_m3 length mismatch"))
        length(E_kg_m3) == cell_count || throw(ArgumentError("E_kg_m3 length mismatch"))
        length(P_mol_m3) == cell_count || throw(ArgumentError("P_mol_m3 length mismatch"))
        clay_system in VALID_CLAY_SYSTEMS || throw(ArgumentError("clay_system must be in VALID_CLAY_SYSTEMS"))

        all(x -> isfinite(x) && x > 0.0, T_K) || throw(ArgumentError("T_K must be finite and > 0 K"))
        all(x -> isfinite(x) && 0.0 < x < 1.0, phi) || throw(ArgumentError("phi must lie strictly in (0, 1)"))
        all(x -> isfinite(x) && x > 0.0, k_m2) || throw(ArgumentError("k_m2 must be strictly positive"))
        all(x -> isfinite(x) && x >= 0.0, B_kg_m3) || throw(ArgumentError("B_kg_m3 must be non-negative"))
        all(x -> isfinite(x) && x >= 0.0, E_kg_m3) || throw(ArgumentError("E_kg_m3 must be non-negative"))
        all(x -> isfinite(x) && x >= 0.0, P_mol_m3) || throw(ArgumentError("P_mol_m3 must be non-negative"))

        for (k, v) in C_aq_mol_m3
            length(v) == cell_count || throw(ArgumentError("C_aq_mol_m3[$k] length mismatch"))
            all(x -> isfinite(x) && x >= 0.0, v) || throw(ArgumentError("C_aq_mol_m3[$k] must be finite and >= 0"))
        end
        for (k, v) in M_min_mol_m3
            length(v) == cell_count || throw(ArgumentError("M_min_mol_m3[$k] length mismatch"))
            all(x -> isfinite(x) && x >= 0.0, v) || throw(ArgumentError("M_min_mol_m3[$k] must be finite and >= 0"))
        end
        for (k, v) in N_iso_mol_m3
            length(v) == cell_count || throw(ArgumentError("N_iso_mol_m3[$k] length mismatch"))
            all(x -> isfinite(x) && x >= 0.0, v) || throw(ArgumentError("N_iso_mol_m3[$k] must be finite and >= 0"))
        end

        new(
            cell_count,
            Float64.(u_m_s),
            Float64.(p_Pa),
            Float64.(T_K),
            Float64.(phi),
            Float64.(k_m2),
            deepcopy(C_aq_mol_m3),
            deepcopy(M_min_mol_m3),
            deepcopy(N_iso_mol_m3),
            Float64.(B_kg_m3),
            Float64.(E_kg_m3),
            Float64.(P_mol_m3),
            clay_system
        )
    end
end

"""
    evaluate_global_mass_conservation_rate(;
        initial_total_mol::Real,
        final_total_mol::Real,
        dt_s::Real,
        J_in_mol_s::Real,
        J_out_mol_s::Real,
        R_nuclear_mol_s::Real = 0.0,
        rel_tol::Real = 1.0e-6)

Verifies the global elemental conservation law:
`dM_i^total / dt = J_i,in - J_i,out + R_i^nuclear`
integrated over `dt_s`.
"""
function evaluate_global_mass_conservation_rate(;
        initial_total_mol::Real,
        final_total_mol::Real,
        dt_s::Real,
        J_in_mol_s::Real,
        J_out_mol_s::Real,
        R_nuclear_mol_s::Real = 0.0,
        rel_tol::Real = 1.0e-6)
    dt_s > 0.0 || throw(ArgumentError("dt_s must be strictly positive"))
    initial_total_mol >= 0.0 || throw(ArgumentError("initial_total_mol must be non-negative"))
    final_total_mol >= 0.0 || throw(ArgumentError("final_total_mol must be non-negative"))

    observed_dM_dt = (Float64(final_total_mol) - Float64(initial_total_mol)) / Float64(dt_s)
    expected_dM_dt = Float64(J_in_mol_s) - Float64(J_out_mol_s) + Float64(R_nuclear_mol_s)
    expected_final_mol = Float64(initial_total_mol) + expected_dM_dt * Float64(dt_s)
    abs_residual_mol = abs(Float64(final_total_mol) - expected_final_mol)
    denom = max(1.0e-15, max(Float64(initial_total_mol), abs(expected_final_mol)))
    rel_residual = abs_residual_mol / denom

    return (
        is_conserved = (rel_residual <= Float64(rel_tol)),
        observed_dM_dt_mol_s = observed_dM_dt,
        expected_dM_dt_mol_s = expected_dM_dt,
        expected_final_mol = expected_final_mol,
        abs_residual_mol = abs_residual_mol,
        rel_residual = rel_residual,
        well_below_one_percent = (rel_residual < 0.01)
    )
end

"""
    evaluate_clay_system_transition(system_id::Symbol; eps_fraction::Real = 0.15, ionic_strength_mol_L::Real = 0.05)

Evaluates the mineralogical and hydraulic transition across the three M7/SMBR clay reference systems:
- `:system_a_sudoite_illite` (Hydrothermal reactor chlorite/sudoite + illite, low swelling, high thermal stability)
- `:system_b_bentonite_montmorillonite` (Swelling smectite/montmorillonite barrier, high CEC, interlayer osmotic swelling)
- `:system_c_illite_smectite_eps_biofilm` (Mixed-layer illite/smectite + EPS/biofilm pore-bridging matrix)
"""
function evaluate_clay_system_transition(
        system_id::Symbol;
        eps_fraction::Real = 0.15,
        ionic_strength_mol_L::Real = 0.05)
    system_id in VALID_CLAY_SYSTEMS || throw(ArgumentError("Unknown clay system :$(system_id)"))
    (0.0 <= eps_fraction < 0.95) || throw(ArgumentError("eps_fraction must lie in [0, 0.95)"))
    ionic_strength_mol_L > 0.0 || throw(ArgumentError("ionic_strength_mol_L must be > 0"))

    if system_id == :system_a_sudoite_illite
        cec_meq_100g = 18.0
        swelling_index = 0.05
        baseline_phi = 0.34
        baseline_k_m2 = 4.2e-13
    elseif system_id == :system_b_bentonite_montmorillonite
        cec_meq_100g = 95.0
        swelling_index = 0.65 / (1.0 + 2.0 * Float64(ionic_strength_mol_L))
        baseline_phi = max(0.12, 0.42 * (1.0 - 0.5 * swelling_index))
        baseline_k_m2 = 1.5e-15 * ((baseline_phi / 0.30) ^ 3)
    else
        cec_meq_100g = 52.0
        swelling_index = 0.28
        baseline_phi = max(0.10, 0.38 * (1.0 - Float64(eps_fraction)))
        baseline_k_m2 = 2.5e-13 * ((baseline_phi / 0.38) ^ 3) * ((0.62 / (1.0 - baseline_phi)) ^ 2)
    end

    return (
        system_id = system_id,
        epistemic_status = :modeled_system,
        cec_meq_100g = cec_meq_100g,
        swelling_index = swelling_index,
        effective_porosity = baseline_phi,
        effective_permeability_m2 = baseline_k_m2
    )
end

const VALID_REACTOR_POOL_IDS = (
    :pool_1_inlet_preconditioning,
    :pool_2_clay_eps_sorption,
    :pool_3_radiolytic_melanin_shield,
    :pool_4_effluent_polishing,
)

"""
    MultiPoolReactorBiofilmState

Cross-session state contract for coupled multi-pool bioreactor hydraulics, radiolysis, and
biofilm/EPS/streamer inventories exchanged between `SmallModularBioReactors.jl`, `M7Multiscale`,
and `MATLAB R2026b` (`Simulink` / `Simscape` 3D EDA).
"""
struct MultiPoolReactorBiofilmState
    pool_ids::Vector{Symbol}
    pool_volumes_L::Vector{Float64}
    inter_pool_flow_L_s::Vector{Float64}
    pool_biomass_attached_kg_m3::Vector{Float64}
    pool_biomass_suspended_kg_m3::Vector{Float64}
    pool_eps_fraction::Vector{Float64}
    pool_porosity_phi::Vector{Float64}
    pool_dose_rate_Gy_s::Vector{Float64}

    function MultiPoolReactorBiofilmState(;
            pool_ids::AbstractVector{Symbol} = collect(VALID_REACTOR_POOL_IDS),
            pool_volumes_L::AbstractVector{<:Real} = [0.25, 0.25, 0.25, 0.25],
            inter_pool_flow_L_s::AbstractVector{<:Real} = [2.5e-4, 2.5e-4, 2.5e-4],
            pool_biomass_attached_kg_m3::AbstractVector{<:Real} = [0.82, 0.74, 0.65, 0.58],
            pool_biomass_suspended_kg_m3::AbstractVector{<:Real} = [0.12, 0.14, 0.11, 0.09],
            pool_eps_fraction::AbstractVector{<:Real} = [0.22, 0.20, 0.18, 0.16],
            pool_porosity_phi::AbstractVector{<:Real} = [0.38, 0.36, 0.35, 0.37],
            pool_dose_rate_Gy_s::AbstractVector{<:Real} = [1.2e-6, 8.5e-7, 4.2e-7, 2.1e-7])
        n = length(pool_ids)
        n >= 2 || throw(ArgumentError("MultiPoolReactorBiofilmState requires at least 2 reactor pools"))
        (length(pool_volumes_L) == n &&
         length(inter_pool_flow_L_s) == n - 1 &&
         length(pool_biomass_attached_kg_m3) == n &&
         length(pool_biomass_suspended_kg_m3) == n &&
         length(pool_eps_fraction) == n &&
         length(pool_porosity_phi) == n &&
         length(pool_dose_rate_Gy_s) == n) ||
            throw(ArgumentError("Dimension mismatch in MultiPoolReactorBiofilmState vectors"))
        all(v -> isfinite(v) && v > 0.0, pool_volumes_L) || throw(ArgumentError("pool_volumes_L must be > 0"))
        all(q -> isfinite(q) && q >= 0.0, inter_pool_flow_L_s) || throw(ArgumentError("inter_pool_flow_L_s must be >= 0"))
        all(b -> isfinite(b) && b >= 0.0, pool_biomass_attached_kg_m3) || throw(ArgumentError("pool_biomass_attached_kg_m3 must be >= 0"))
        all(b -> isfinite(b) && b >= 0.0, pool_biomass_suspended_kg_m3) || throw(ArgumentError("pool_biomass_suspended_kg_m3 must be >= 0"))
        all(e -> isfinite(e) && 0.0 <= e < 0.95, pool_eps_fraction) || throw(ArgumentError("pool_eps_fraction must lie in [0, 0.95)"))
        all(p -> isfinite(p) && 0.0 < p < 1.0, pool_porosity_phi) || throw(ArgumentError("pool_porosity_phi must lie in (0, 1)"))
        all(d -> isfinite(d) && d >= 0.0, pool_dose_rate_Gy_s) || throw(ArgumentError("pool_dose_rate_Gy_s must be >= 0"))

        new(
            Symbol.(pool_ids),
            Float64.(pool_volumes_L),
            Float64.(inter_pool_flow_L_s),
            Float64.(pool_biomass_attached_kg_m3),
            Float64.(pool_biomass_suspended_kg_m3),
            Float64.(pool_eps_fraction),
            Float64.(pool_porosity_phi),
            Float64.(pool_dose_rate_Gy_s)
        )
    end
end

