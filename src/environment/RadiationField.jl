# ==============================================================================
# src/environment/RadiationField.jl
# Spectral Dose Rate Integration D_dot(x,t) = int Phi(E,x,t) * E_dep(E,x,t) dE,
# Spatial Dose Gradient grad(D_dot), and Radiolytic Species Generation
# (Strictly Separated from SPICE Ephemeris Geometry)
# ==============================================================================

struct SpectralRadiationFieldSpec
    energy_bins_MeV::Vector{Float64}
    incident_flux_particles_m2_s_MeV::Vector{Float64}
    mass_energy_absorption_m2_kg::Vector{Float64}
    linear_attenuation_m_inv::Vector{Float64}

    function SpectralRadiationFieldSpec(;
            energy_bins_MeV::AbstractVector{<:Real} = [0.5, 1.5, 5.0, 20.0, 100.0],
            incident_flux_particles_m2_s_MeV::AbstractVector{<:Real} = [2.5e7, 1.2e7, 4.0e6, 8.0e5, 1.5e5],
            mass_energy_absorption_m2_kg::AbstractVector{<:Real} = [0.0030, 0.0028, 0.0024, 0.0020, 0.0018],
            linear_attenuation_m_inv::AbstractVector{<:Real} = [18.0, 12.0, 7.5, 4.2, 2.5])
        n = length(energy_bins_MeV)
        n >= 2 || throw(ArgumentError("energy_bins_MeV must have at least 2 bins for quadrature"))
        (length(incident_flux_particles_m2_s_MeV) == n &&
         length(mass_energy_absorption_m2_kg) == n &&
         length(linear_attenuation_m_inv) == n) ||
            throw(ArgumentError("All spectral vectors must have matching length $(n)"))
        issorted(energy_bins_MeV) || throw(ArgumentError("energy_bins_MeV must be sorted in ascending order"))
        all(x -> isfinite(x) && x > 0.0, energy_bins_MeV) || throw(ArgumentError("energy_bins_MeV must be > 0"))
        all(x -> isfinite(x) && x >= 0.0, incident_flux_particles_m2_s_MeV) ||
            throw(ArgumentError("incident_flux_particles_m2_s_MeV must be >= 0"))
        all(x -> isfinite(x) && x > 0.0, mass_energy_absorption_m2_kg) ||
            throw(ArgumentError("mass_energy_absorption_m2_kg must be > 0"))
        all(x -> isfinite(x) && x >= 0.0, linear_attenuation_m_inv) ||
            throw(ArgumentError("linear_attenuation_m_inv must be >= 0"))

        new(
            Float64.(energy_bins_MeV),
            Float64.(incident_flux_particles_m2_s_MeV),
            Float64.(mass_energy_absorption_m2_kg),
            Float64.(linear_attenuation_m_inv)
        )
    end
end

"""
    compute_dose_rate_profile(spec::SpectralRadiationFieldSpec, depths_m::AbstractVector{<:Real};
        melanin_mass_fraction::Real = 0.0,
        shielding_thickness_m::Real = 0.005)

Evaluates the spectral integral:
`D_dot(x, t) = int Phi(E, x, t) * E_dep(E, x, t) dE` (in `Gy/s` and `Gy/hr`)
across cell depths `depths_m`, computing the spatial dose-rate gradient `d(D_dot)/dx` (`Gy/(s*m)`)
and primary water-radiolysis generation rates (`H2`, `H2O2`, `OH_radical` in `mol/(m^3*s)`).
"""
function compute_dose_rate_profile(
        spec::SpectralRadiationFieldSpec,
        depths_m::AbstractVector{<:Real};
        melanin_mass_fraction::Real = 0.0,
        shielding_thickness_m::Real = 0.005,
        fluid_density_kg_m3::Real = 1000.0)
    n_cells = length(depths_m)
    n_cells >= 2 || throw(ArgumentError("depths_m must contain at least 2 spatial points"))
    all(x -> isfinite(x) && x >= 0.0, depths_m) || throw(ArgumentError("depths_m must be finite and >= 0"))
    issorted(depths_m; lt = <=) || throw(ArgumentError("depths_m must be strictly increasing"))
    (0.0 <= melanin_mass_fraction <= 0.5) || throw(ArgumentError("melanin_mass_fraction must lie in [0, 0.5]"))
    shielding_thickness_m >= 0.0 || throw(ArgumentError("shielding_thickness_m must be >= 0"))
    fluid_density_kg_m3 > 0.0 || throw(ArgumentError("fluid_density_kg_m3 must be > 0"))

    mev_to_joule = 1.602176634e-13
    E_MeV = spec.energy_bins_MeV
    n_bins = length(E_MeV)
    # Melanin enhances linear attenuation of ionizing/low-energy secondary photons/electrons
    melanin_att_boost = 1.0 + 1.8 * Float64(melanin_mass_fraction)

    dose_rate_Gy_s = zeros(Float64, n_cells)
    for i in 1:n_cells
        total_depth = Float64(shielding_thickness_m) + Float64(depths_m[i])
        integrand = zeros(Float64, n_bins)
        for k in 1:n_bins
            att = exp(-spec.linear_attenuation_m_inv[k] * melanin_att_boost * total_depth)
            phi_E = spec.incident_flux_particles_m2_s_MeV[k] * att
            # E_dep per unit mass (J / kg = Gy) per unit flux: E_MeV * mev_to_joule * mu_en/rho
            e_dep_J_kg = E_MeV[k] * mev_to_joule * spec.mass_energy_absorption_m2_kg[k]
            integrand[k] = phi_E * e_dep_J_kg
        end
        # Trapezoidal integration over energy bins
        s = 0.0
        for k in 1:(n_bins - 1)
            dE = E_MeV[k + 1] - E_MeV[k]
            s += 0.5 * (integrand[k] + integrand[k + 1]) * dE
        end
        dose_rate_Gy_s[i] = s
    end

    # Finite-difference spatial gradient grad_D_dot (Gy / (s * m))
    grad_dose_rate_Gy_s_m = zeros(Float64, n_cells)
    for i in 1:n_cells
        if i == 1
            dx = max(1.0e-9, Float64(depths_m[2] - depths_m[1]))
            grad_dose_rate_Gy_s_m[i] = (dose_rate_Gy_s[2] - dose_rate_Gy_s[1]) / dx
        elseif i == n_cells
            dx = max(1.0e-9, Float64(depths_m[n_cells] - depths_m[n_cells - 1]))
            grad_dose_rate_Gy_s_m[i] = (dose_rate_Gy_s[n_cells] - dose_rate_Gy_s[n_cells - 1]) / dx
        else
            dx = max(1.0e-9, Float64(depths_m[i + 1] - depths_m[i - 1]))
            grad_dose_rate_Gy_s_m[i] = (dose_rate_Gy_s[i + 1] - dose_rate_Gy_s[i - 1]) / dx
        end
    end

    # Radiolytic G-values (mol / J): H2 ~ 0.45e-7, H2O2 ~ 0.70e-7, OH_radical ~ 2.80e-7
    # Melanin radical scavenging reduces free OH_radical / H2O2 in bulk pore water
    scavenge_factor = 1.0 / (1.0 + 4.5 * Float64(melanin_mass_fraction))
    rho = Float64(fluid_density_kg_m3)
    R_H2_mol_m3_s = dose_rate_Gy_s .* (rho * 0.45e-7)
    R_H2O2_mol_m3_s = dose_rate_Gy_s .* (rho * 0.70e-7 * scavenge_factor)
    R_OH_mol_m3_s = dose_rate_Gy_s .* (rho * 2.80e-7 * scavenge_factor)

    return (
        depths_m = Float64.(depths_m),
        dose_rate_Gy_s = dose_rate_Gy_s,
        dose_rate_Gy_hr = dose_rate_Gy_s .* 3600.0,
        grad_dose_rate_Gy_s_m = grad_dose_rate_Gy_s_m,
        radiolytic_H2_mol_m3_s = R_H2_mol_m3_s,
        radiolytic_H2O2_mol_m3_s = R_H2O2_mol_m3_s,
        radiolytic_OH_mol_m3_s = R_OH_mol_m3_s,
        spice_geometry_decoupled = true
    )
end
