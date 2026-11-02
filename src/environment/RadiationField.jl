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

"""
    HealthPhysicsALARAContract

Health-Physics ALARA (`As Low As Reasonably Achievable`: Time `Delta_t`, Distance `1/r^2`,
Shielding `exp(-mu*z)` + biological melanin/Neurospora radioprotection) and Epistemic/Operational
Hippocratism (`primum non nocere` / first, do no harm) contract for coupled nuclear/radiobiological
and bioreactor systems.
"""
struct HealthPhysicsALARAContract
    time_exposure_s::Float64
    max_allowable_cumulative_dose_Gy::Float64
    reference_distance_m::Float64
    standoff_distance_m::Float64
    shielding_regolith_thickness_m::Float64
    shielding_water_clay_thickness_m::Float64
    melanin_mass_fraction::Float64
    neurospora_circadian_gate_active::Bool
    epistemic_hippocratism_enforced::Bool

    function HealthPhysicsALARAContract(;
            time_exposure_s::Real = 3600.0,
            max_allowable_cumulative_dose_Gy::Real = 0.50,
            reference_distance_m::Real = 1.0,
            standoff_distance_m::Real = 10.0,
            shielding_regolith_thickness_m::Real = 0.08,
            shielding_water_clay_thickness_m::Real = 0.02,
            melanin_mass_fraction::Real = 0.18,
            neurospora_circadian_gate_active::Bool = true,
            epistemic_hippocratism_enforced::Bool = true)
        (isfinite(time_exposure_s) && time_exposure_s >= 0.0) ||
            throw(ArgumentError("time_exposure_s must be finite and >= 0"))
        (isfinite(max_allowable_cumulative_dose_Gy) && max_allowable_cumulative_dose_Gy > 0.0) ||
            throw(ArgumentError("max_allowable_cumulative_dose_Gy must be > 0"))
        (isfinite(reference_distance_m) && reference_distance_m > 0.0) ||
            throw(ArgumentError("reference_distance_m must be > 0"))
        (isfinite(standoff_distance_m) && standoff_distance_m >= reference_distance_m) ||
            throw(ArgumentError("standoff_distance_m must be >= reference_distance_m for 1/r^2 ALARA separation"))
        (isfinite(shielding_regolith_thickness_m) && shielding_regolith_thickness_m >= 0.0) ||
            throw(ArgumentError("shielding_regolith_thickness_m must be >= 0"))
        (isfinite(shielding_water_clay_thickness_m) && shielding_water_clay_thickness_m >= 0.0) ||
            throw(ArgumentError("shielding_water_clay_thickness_m must be >= 0"))
        (isfinite(melanin_mass_fraction) && 0.0 <= melanin_mass_fraction <= 0.50) ||
            throw(ArgumentError("melanin_mass_fraction must lie in [0.0, 0.50]"))
        epistemic_hippocratism_enforced ||
            throw(ArgumentError("Health-Physics Hippocratism Violation: epistemic_hippocratism_enforced must be true (primum non nocere)"))

        new(
            Float64(time_exposure_s),
            Float64(max_allowable_cumulative_dose_Gy),
            Float64(reference_distance_m),
            Float64(standoff_distance_m),
            Float64(shielding_regolith_thickness_m),
            Float64(shielding_water_clay_thickness_m),
            Float64(melanin_mass_fraction),
            neurospora_circadian_gate_active,
            epistemic_hippocratism_enforced
        )
    end
end

"""
    evaluate_health_physics_alara_and_hippocratism(
        spec::SpectralRadiationFieldSpec = SpectralRadiationFieldSpec(),
        contract::HealthPhysicsALARAContract = HealthPhysicsALARAContract();
        depths_m::AbstractVector{<:Real} = [0.0, 0.002, 0.004, 0.006],
        temperature_K::Real = 298.15,
        spice_units_verified_by_spkpos::Bool = true,
        claim_unverified_low_dose_hormesis_as_safe::Bool = false,
        recursion_depth::Integer = 0)

Evaluates the 3-pillar Health-Physics ALARA state (`Time Delta_t`, `Distance 1/r^2`,
`Shielding exp(-mu*z)` + `Melanin/Neurospora` biological radical scavenging) together with
Epistemic & Operational Hippocratism (`primum non nocere` / first, do no harm):
- Rejects unestablished low-dose hormesis / net radiosynthesis claims used to relax shielding (`ArgumentError`).
- Rejects recursive self-modification (`recursion_depth != 0` throws `ArgumentError`).
- Holds out-of-range thermal regimes (`T < 280 K`) or unverified `SPICE.jl` unit assumptions in
  `:Incomplete` status rather than silently promoting them to `:Accepted`.
"""
function evaluate_health_physics_alara_and_hippocratism(
        spec::SpectralRadiationFieldSpec = SpectralRadiationFieldSpec(),
        contract::HealthPhysicsALARAContract = HealthPhysicsALARAContract();
        depths_m::AbstractVector{<:Real} = [0.0, 0.002, 0.004, 0.006],
        temperature_K::Real = 298.15,
        spice_units_verified_by_spkpos::Bool = true,
        claim_unverified_low_dose_hormesis_as_safe::Bool = false,
        recursion_depth::Integer = 0)
    claim_unverified_low_dose_hormesis_as_safe && throw(ArgumentError(
        "Health-Physics Hippocratism Violation (Primum Non Nocere): unestablished low-dose radiation hormesis or direct ATP radiosynthesis cannot be used to relax ALARA shielding or exposure limits"
    ))
    recursion_depth == 0 || throw(ArgumentError(
        "Operational Hippocratism Violation: non-recursive consolidation requires recursion_depth == 0, got $(recursion_depth)"
    ))
    (isfinite(temperature_K) && temperature_K > 0.0) ||
        throw(ArgumentError("temperature_K must be finite and > 0"))

    # 1. Baseline unshielded reference dose rate at reference_distance_m
    unshielded_ref = compute_dose_rate_profile(
        spec,
        depths_m;
        melanin_mass_fraction = 0.0,
        shielding_thickness_m = 0.0
    )

    # 2. Distance (1/r^2) + Shielding (regolith + water/clay + melanin) attenuation
    inv_sq_factor = (contract.reference_distance_m^2) / (contract.standoff_distance_m^2)
    total_shielding_m = contract.shielding_regolith_thickness_m + contract.shielding_water_clay_thickness_m
    shielded_raw = compute_dose_rate_profile(
        spec,
        depths_m;
        melanin_mass_fraction = contract.melanin_mass_fraction,
        shielding_thickness_m = total_shielding_m
    )

    alara_dose_rate_Gy_s = shielded_raw.dose_rate_Gy_s .* inv_sq_factor
    alara_dose_rate_Gy_hr = alara_dose_rate_Gy_s .* 3600.0
    cumulative_dose_Gy = alara_dose_rate_Gy_s .* contract.time_exposure_s

    peak_unshielded_Gy_s = maximum(unshielded_ref.dose_rate_Gy_s)
    peak_alara_Gy_s = maximum(alara_dose_rate_Gy_s)
    peak_cumulative_Gy = maximum(cumulative_dose_Gy)
    max_safe_exposure_time_s = contract.max_allowable_cumulative_dose_Gy / max(peak_alara_Gy_s, 1.0e-30)
    alara_dose_reduction_factor = peak_unshielded_Gy_s / max(peak_alara_Gy_s, 1.0e-30)

    # Automated payload safing triggers at 80% of max_allowable_cumulative_dose_Gy (conservative ALARA margin)
    payload_safing_triggered = peak_cumulative_Gy >= 0.80 * contract.max_allowable_cumulative_dose_Gy
    within_alara_dose_budget = peak_cumulative_Gy <= contract.max_allowable_cumulative_dose_Gy

    # 3. Epistemic Hippocratism 3-valued gate
    below_280K_underconstrained = Float64(temperature_K) < 280.0
    epistemic_verdict = if !within_alara_dose_budget
        :Rejected
    elseif below_280K_underconstrained || !spice_units_verified_by_spkpos || payload_safing_triggered
        :Incomplete
    else
        :Accepted
    end

    return (
        alara_satisfied = within_alara_dose_budget && (alara_dose_reduction_factor > 10.0),
        epistemic_hippocratism_verdict = epistemic_verdict,
        primum_non_nocere_enforced = contract.epistemic_hippocratism_enforced && (recursion_depth == 0),
        time_control = (
            exposure_time_s = contract.time_exposure_s,
            max_safe_exposure_time_s = max_safe_exposure_time_s,
            peak_cumulative_dose_Gy = peak_cumulative_Gy,
            max_allowable_cumulative_dose_Gy = contract.max_allowable_cumulative_dose_Gy,
            payload_safing_triggered = payload_safing_triggered
        ),
        distance_control = (
            reference_distance_m = contract.reference_distance_m,
            standoff_distance_m = contract.standoff_distance_m,
            inverse_square_attenuation_factor = inv_sq_factor
        ),
        shielding_control = (
            regolith_thickness_m = contract.shielding_regolith_thickness_m,
            water_clay_thickness_m = contract.shielding_water_clay_thickness_m,
            total_physical_shielding_m = total_shielding_m,
            melanin_mass_fraction = contract.melanin_mass_fraction,
            alara_dose_reduction_factor = alara_dose_reduction_factor
        ),
        dose_rate_Gy_s = alara_dose_rate_Gy_s,
        dose_rate_Gy_hr = alara_dose_rate_Gy_hr,
        cumulative_dose_Gy = cumulative_dose_Gy,
        below_280K_underconstrained = below_280K_underconstrained,
        spice_units_verified_by_spkpos = spice_units_verified_by_spkpos,
        recursion_depth = Int(recursion_depth)
    )
end

