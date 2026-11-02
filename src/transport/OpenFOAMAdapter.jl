# ==============================================================================
# src/transport/OpenFOAMAdapter.jl
# Mass-Conservative Darcy-Brinkman-Stokes (DBS) Porous Media Flow & Reactive
# Transport Adapter (Menke et al. 2026 GeoChemFoam & arXiv:2610.01565)
# ==============================================================================

"""
    compute_kozeny_carman_permeability(phi::Real; phi0::Real = 0.38, k0_m2::Real = 2.5e-13) -> Float64

Computes the local porous-medium permeability `k(phi)` (`m^2`) as biofilm/EPS or mineral
precipitation modifies porosity `phi`.
"""
function compute_kozeny_carman_permeability(phi::Real; phi0::Real = 0.38, k0_m2::Real = 2.5e-13)
    (0.0 < phi < 1.0 && 0.0 < phi0 < 1.0 && k0_m2 > 0.0) ||
        throw(ArgumentError("Porosities phi, phi0 must lie in (0, 1) and k0_m2 must be > 0"))
    p = Float64(phi)
    p0 = Float64(phi0)
    return Float64(k0_m2) * ((p / p0) ^ 3) * (((1.0 - p0) / (1.0 - p)) ^ 2)
end

"""
    step_darcy_brinkman_stokes_transport(;
        C_aq_mol_m3::AbstractVector{<:Real},
        phi::AbstractVector{<:Real},
        u_darcy_m_s::Real = 2.0e-5,
        diffusivity_m2_s::Real = 1.0e-9,
        dx_m::Real = 1.0e-3,
        dt_s::Real = 5.0,
        closed_boundaries::Bool = true)

Executes a mass-conservative finite-volume advection-diffusion step for dissolved pore-water
species in a Darcy-Brinkman-Stokes (DBS) porous medium (`Menke et al. 2026`).
When `closed_boundaries == true` (zero-flux Neumann ends), total pore-water moles
`sum(phi_i * C_aq_i * V_i)` are conserved to machine precision (`residual << 1%`).
"""
function step_darcy_brinkman_stokes_transport(;
        C_aq_mol_m3::AbstractVector{<:Real},
        phi::AbstractVector{<:Real},
        u_darcy_m_s::Real = 2.0e-5,
        diffusivity_m2_s::Real = 1.0e-9,
        dx_m::Real = 1.0e-3,
        dt_s::Real = 5.0,
        cell_area_m2::Real = 1.0e-6,
        closed_boundaries::Bool = true)
    n = length(C_aq_mol_m3)
    n >= 2 || throw(ArgumentError("C_aq_mol_m3 must contain at least 2 cells"))
    length(phi) == n || throw(ArgumentError("phi length must match C_aq_mol_m3"))
    all(c -> isfinite(c) && c >= 0.0, C_aq_mol_m3) || throw(ArgumentError("C_aq_mol_m3 must be finite and >= 0"))
    all(p -> isfinite(p) && 0.0 < p < 1.0, phi) || throw(ArgumentError("phi must lie in (0, 1)"))
    (diffusivity_m2_s > 0.0 && dx_m > 0.0 && dt_s > 0.0 && cell_area_m2 > 0.0) ||
        throw(ArgumentError("diffusivity_m2_s, dx_m, dt_s, and cell_area_m2 must be > 0"))

    u = Float64(u_darcy_m_s)
    D = Float64(diffusivity_m2_s)
    dx = Float64(dx_m)
    dt = Float64(dt_s)
    vol = dx * Float64(cell_area_m2)

    cfl = (abs(u) * dt) / (minimum(phi) * dx)
    fourier = (D * dt) / (minimum(phi) * (dx ^ 2))
    peclet = (abs(u) * dx) / D
    (cfl + 2.0 * fourier <= 0.95) ||
        throw(ArgumentError("Explicit DBS transport step violates CFL/Fourier stability (cfl=$(cfl), fourier=$(fourier))"))

    # Compute conservative face fluxes F_{f} (mol / (m^2 * s)) across interior faces f = 1..(n-1)
    face_flux = zeros(Float64, n + 1)
    if !closed_boundaries
        face_flux[1] = u * Float64(C_aq_mol_m3[1])
        face_flux[n + 1] = u * Float64(C_aq_mol_m3[n])
    end

    for f in 2:n
        iL = f - 1
        iR = f
        cL = Float64(C_aq_mol_m3[iL])
        cR = Float64(C_aq_mol_m3[iR])
        phi_f = 0.5 * (Float64(phi[iL]) + Float64(phi[iR]))
        adv = u >= 0.0 ? u * cL : u * cR
        diff = -phi_f * D * (cR - cL) / dx
        face_flux[f] = adv + diff
    end

    C_next = zeros(Float64, n)
    for i in 1:n
        net_flux_div = (face_flux[i] - face_flux[i + 1]) / dx
        moles_per_bulk_vol = Float64(phi[i]) * Float64(C_aq_mol_m3[i]) + dt * net_flux_div
        C_next[i] = max(0.0, moles_per_bulk_vol / Float64(phi[i]))
    end

    initial_moles = sum(Float64(phi[i]) * Float64(C_aq_mol_m3[i]) * vol for i in 1:n)
    final_moles = sum(Float64(phi[i]) * C_next[i] * vol for i in 1:n)
    boundary_net_moles = (face_flux[1] - face_flux[n + 1]) * Float64(cell_area_m2) * dt
    expected_final_moles = initial_moles + boundary_net_moles
    rel_mass_residual = abs(final_moles - expected_final_moles) / max(1.0e-18, initial_moles)

    return (
        C_next_mol_m3 = C_next,
        cfl_number = cfl,
        peclet_number = peclet,
        initial_moles = initial_moles,
        final_moles = final_moles,
        rel_mass_residual = rel_mass_residual,
        mass_conserved_below_1pct = (rel_mass_residual < 0.01)
    )
end
