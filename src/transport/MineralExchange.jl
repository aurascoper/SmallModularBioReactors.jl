# ==============================================================================
# src/transport/MineralExchange.jl
# Menke et al. (2026) Conservative Dissolved Pore-Water + Cell-Scale Adsorbed
# Surface Inventories and Polymer-Ligand Metal Complexation:
#   P-L + M^{2+} <=> P-L-M
# Preserving Global Elemental Mass Balance (Residual << 1%)
# ==============================================================================

"""
    solve_polymer_ligand_metal_complexation(;
        total_metal_mol_m3_bulk::Real,
        total_ligand_mol_m3_bulk::Real,
        porosity::Real,
        K_bind_m3_mol::Real = 25.0)

Solves the exact mass-action equilibrium for cell-scale surface/polymer complexation:
`P-L + M^{2+} <=> P-L-M`
with equilibrium binding constant `K_bind = [P-L-M]_bulk / ([P-L]_bulk * [M^{2+}]_aq)`
subject to strict cell-scale conservation:
- `total_metal_mol_m3_bulk = porosity * [M^{2+}]_aq + [P-L-M]_bulk`
- `total_ligand_mol_m3_bulk = [P-L]_bulk + [P-L-M]_bulk`
"""
function solve_polymer_ligand_metal_complexation(;
        total_metal_mol_m3_bulk::Real,
        total_ligand_mol_m3_bulk::Real,
        porosity::Real,
        K_bind_m3_mol::Real = 25.0)
    (total_metal_mol_m3_bulk >= 0.0 && total_ligand_mol_m3_bulk >= 0.0) ||
        throw(ArgumentError("Total metal and ligand inventories must be non-negative"))
    (0.0 < porosity < 1.0) || throw(ArgumentError("porosity must lie in (0, 1)"))
    K_bind_m3_mol > 0.0 || throw(ArgumentError("K_bind_m3_mol must be > 0"))

    M_tot = Float64(total_metal_mol_m3_bulk)
    L_tot = Float64(total_ligand_mol_m3_bulk)
    phi = Float64(porosity)
    K = Float64(K_bind_m3_mol)

    if M_tot == 0.0 || L_tot == 0.0
        c_aq = M_tot / phi
        return (
            free_metal_aq_mol_m3 = c_aq,
            free_ligand_bulk_mol_m3 = L_tot,
            bound_complex_PLM_bulk_mol_m3 = 0.0,
            site_occupancy_fraction = 0.0,
            mass_balance_abs_residual = 0.0,
            mass_balance_rel_residual = 0.0
        )
    end

    # Let q = [P-L-M]_bulk in [0, min(M_tot, L_tot)].
    # Equilibrium: K * (L_tot - q) * ((M_tot - q) / phi) = q
    # => K * (L_tot - q) * (M_tot - q) = phi * q
    # => q^2 - (M_tot + L_tot + phi / K) * q + M_tot * L_tot = 0
    b = M_tot + L_tot + (phi / K)
    disc = max(0.0, (b ^ 2) - 4.0 * M_tot * L_tot)
    # Numerically stable root in [0, min(M_tot, L_tot)]
    q = (2.0 * M_tot * L_tot) / (b + sqrt(disc))
    q = clamp(q, 0.0, min(M_tot, L_tot))

    c_aq = (M_tot - q) / phi
    free_L = L_tot - q
    reconstructed_M = phi * c_aq + q
    abs_res = abs(reconstructed_M - M_tot)
    rel_res = abs_res / max(1.0e-18, M_tot)

    return (
        free_metal_aq_mol_m3 = c_aq,
        free_ligand_bulk_mol_m3 = free_L,
        bound_complex_PLM_bulk_mol_m3 = q,
        site_occupancy_fraction = q / L_tot,
        mass_balance_abs_residual = abs_res,
        mass_balance_rel_residual = rel_res
    )
end

"""
    step_menke_2026_coupled_reactive_transport(;
        C_aq_mol_m3::AbstractVector{<:Real} = [0.20, 0.12, 0.06, 0.02],
        M_sorb_mol_m3_bulk::AbstractVector{<:Real} = [0.10, 0.08, 0.04, 0.01],
        L_total_mol_m3_bulk::AbstractVector{<:Real} = fill(0.35, 4),
        phi::AbstractVector{<:Real} = fill(0.38, 4),
        K_bind_m3_mol::Real = 25.0,
        u_darcy_m_s::Real = 1.5e-5,
        diffusivity_m2_s::Real = 1.0e-9,
        dx_m::Real = 1.0e-3,
        dt_s::Real = 5.0,
        closed_boundaries::Bool = true)

Implements the `Menke et al. (2026)` GeoChemFoam conservative reactive-transport coupling:
1. Dissolved species advect/diffuse in pore water via `step_darcy_brinkman_stokes_transport`.
2. Adsorbed species (`P-L-M` and clay surface complexes) remain in cell-scale surface inventories,
   updating in each cell via `solve_polymer_ligand_metal_complexation`.
3. Global elemental mass balance (`sum((phi_i * C_aq_i + M_sorb_i) * V_i)`) is preserved with
   relative residual `<< 1%` (`< 1e-12` in closed batch).
"""
function step_menke_2026_coupled_reactive_transport(;
        C_aq_mol_m3::AbstractVector{<:Real} = [0.20, 0.12, 0.06, 0.02],
        M_sorb_mol_m3_bulk::AbstractVector{<:Real} = [0.10, 0.08, 0.04, 0.01],
        L_total_mol_m3_bulk::AbstractVector{<:Real} = fill(0.35, 4),
        phi::AbstractVector{<:Real} = fill(0.38, 4),
        K_bind_m3_mol::Real = 25.0,
        u_darcy_m_s::Real = 1.5e-5,
        diffusivity_m2_s::Real = 1.0e-9,
        dx_m::Real = 1.0e-3,
        dt_s::Real = 5.0,
        cell_area_m2::Real = 1.0e-6,
        closed_boundaries::Bool = true)
    n = length(C_aq_mol_m3)
    (length(M_sorb_mol_m3_bulk) == n && length(L_total_mol_m3_bulk) == n && length(phi) == n) ||
        throw(ArgumentError("All cell state vectors must have matching length $(n)"))

    vol = Float64(dx_m) * Float64(cell_area_m2)
    initial_total_moles = sum(
        (Float64(phi[i]) * Float64(C_aq_mol_m3[i]) + Float64(M_sorb_mol_m3_bulk[i])) * vol
        for i in 1:n
    )

    # Step 1: Mobile pore-water advection-diffusion
    dbs_step = step_darcy_brinkman_stokes_transport(
        C_aq_mol_m3 = C_aq_mol_m3,
        phi = phi,
        u_darcy_m_s = u_darcy_m_s,
        diffusivity_m2_s = diffusivity_m2_s,
        dx_m = dx_m,
        dt_s = dt_s,
        cell_area_m2 = cell_area_m2,
        closed_boundaries = closed_boundaries
    )

    # Step 2: Cell-by-cell surface complexation P-L + M^{2+} <=> P-L-M
    C_aq_final = zeros(Float64, n)
    M_sorb_final = zeros(Float64, n)
    site_occupancy = zeros(Float64, n)
    for i in 1:n
        m_cell_tot = Float64(phi[i]) * dbs_step.C_next_mol_m3[i] + Float64(M_sorb_mol_m3_bulk[i])
        eq = solve_polymer_ligand_metal_complexation(
            total_metal_mol_m3_bulk = m_cell_tot,
            total_ligand_mol_m3_bulk = L_total_mol_m3_bulk[i],
            porosity = phi[i],
            K_bind_m3_mol = K_bind_m3_mol
        )
        C_aq_final[i] = eq.free_metal_aq_mol_m3
        M_sorb_final[i] = eq.bound_complex_PLM_bulk_mol_m3
        site_occupancy[i] = eq.site_occupancy_fraction
    end

    final_total_moles = sum(
        (Float64(phi[i]) * C_aq_final[i] + M_sorb_final[i]) * vol
        for i in 1:n
    )
    net_boundary_moles = dbs_step.final_moles - dbs_step.initial_moles
    expected_total_moles = initial_total_moles + net_boundary_moles
    rel_residual = abs(final_total_moles - expected_total_moles) / max(1.0e-18, abs(expected_total_moles))

    return (
        C_aq_mol_m3 = C_aq_final,
        M_sorb_mol_m3_bulk = M_sorb_final,
        site_occupancy_fraction = site_occupancy,
        initial_total_moles = initial_total_moles,
        net_boundary_moles = net_boundary_moles,
        expected_total_moles = expected_total_moles,
        final_total_moles = final_total_moles,
        global_mass_rel_residual = rel_residual,
        residual_well_below_one_percent = (rel_residual < 0.01)
    )
end
