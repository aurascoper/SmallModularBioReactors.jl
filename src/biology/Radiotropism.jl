# ==============================================================================
# src/biology/Radiotropism.jl
# Radiation Exposure Response (mu_growth, mu_damage, mu_lysis) vs Candidate
# Directional Radiotropism Flux J_B = -D_B * grad(B) + chi_R * B * grad(D_dot)
# Gated by Mandatory 4-Control Artifact Battery
# ==============================================================================

struct RadiotropismControlBattery
    differential_killing_controlled::Bool
    thermal_gradient_controlled::Bool
    chemical_gradient_controlled::Bool
    optical_mechanical_artifact_controlled::Bool
    confounder_refutation_detected::Bool
end

"""
    evaluate_radiotropism_flux_and_controls(;
        B_kg_m3::AbstractVector{<:Real} = [0.40, 0.45, 0.50, 0.55],
        dose_rate_Gy_s::AbstractVector{<:Real} = [2.0e-4, 1.4e-4, 0.9e-4, 0.5e-4],
        grad_dose_rate_Gy_s_m::AbstractVector{<:Real} = fill(-0.05, 4),
        dx_m::Real = 1.0e-3,
        D_B_m2_s::Real = 2.0e-11,
        chi_R_m2_Gy::Real = 1.5e-6,
        controls::RadiotropismControlBattery = RadiotropismControlBattery(true, true, true, true, false))

Evaluates:
1. Mechanism A (Scalar dose-rate response): `mu_growth(D_dot)`, `mu_damage(D_dot)`, `mu_lysis(D_dot)`.
2. Mechanism B (Directional radiotropism candidate flux): `J_B = -D_B * grad(B) + chi_R * B * grad(D_dot)`.
3. Mandatory 4-control battery (differential killing, thermal gradients, chemical gradients, optical/mechanical artifacts).
   Directional radiotropism is only `:Accepted` if all 4 controls pass and `confounder_refutation_detected == false`.
   If controls are missing without refutation, verdict is `:Incomplete` (never collapsed into `:Rejected`).
"""
function evaluate_radiotropism_flux_and_controls(;
        B_kg_m3::AbstractVector{<:Real} = [0.40, 0.45, 0.50, 0.55],
        dose_rate_Gy_s::AbstractVector{<:Real} = [2.0e-4, 1.4e-4, 0.9e-4, 0.5e-4],
        grad_dose_rate_Gy_s_m::AbstractVector{<:Real} = fill(-0.05, 4),
        dx_m::Real = 1.0e-3,
        D_B_m2_s::Real = 2.0e-11,
        chi_R_m2_Gy::Real = 1.5e-6,
        controls::RadiotropismControlBattery = RadiotropismControlBattery(true, true, true, true, false))
    n = length(B_kg_m3)
    n >= 2 || throw(ArgumentError("B_kg_m3 must have at least 2 cells"))
    (length(dose_rate_Gy_s) == n && length(grad_dose_rate_Gy_s_m) == n) ||
        throw(ArgumentError("B_kg_m3, dose_rate_Gy_s, and grad_dose_rate_Gy_s_m must have matching length"))
    (dx_m > 0.0 && D_B_m2_s > 0.0) || throw(ArgumentError("dx_m and D_B_m2_s must be > 0"))

    dx = Float64(dx_m)
    D_B = Float64(D_B_m2_s)
    chi_R = Float64(chi_R_m2_Gy)

    # Mechanism A: Scalar dose-rate response functions (1/s)
    mu_growth_s_inv = [2.5e-5 / (1.0 + 500.0 * Float64(d)) for d in dose_rate_Gy_s]
    mu_damage_s_inv = [1.2e-2 * Float64(d) for d in dose_rate_Gy_s]
    mu_lysis_s_inv = [4.0e-3 * (Float64(d) ^ 1.2) for d in dose_rate_Gy_s]

    # Mechanism B: Directional flux J_B = -D_B * grad(B) + chi_R * B * grad(D_dot)
    grad_B = zeros(Float64, n)
    for i in 1:n
        if i == 1
            grad_B[i] = (Float64(B_kg_m3[2]) - Float64(B_kg_m3[1])) / dx
        elseif i == n
            grad_B[i] = (Float64(B_kg_m3[n]) - Float64(B_kg_m3[n - 1])) / dx
        else
            grad_B[i] = (Float64(B_kg_m3[i + 1]) - Float64(B_kg_m3[i - 1])) / (2.0 * dx)
        end
    end

    diffusive_flux = -D_B .* grad_B
    radiotropic_drift_flux = [chi_R * Float64(B_kg_m3[i]) * Float64(grad_dose_rate_Gy_s_m[i]) for i in 1:n]
    J_B_total = diffusive_flux .+ radiotropic_drift_flux

    all_controls_passed = controls.differential_killing_controlled &&
                          controls.thermal_gradient_controlled &&
                          controls.chemical_gradient_controlled &&
                          controls.optical_mechanical_artifact_controlled

    verdict = if controls.confounder_refutation_detected
        :Rejected
    elseif all_controls_passed
        :Accepted
    else
        :Incomplete
    end

    return (
        verdict = verdict,
        all_four_controls_passed = all_controls_passed,
        controls = controls,
        scalar_dose_response = (
            mu_growth_s_inv = mu_growth_s_inv,
            mu_damage_s_inv = mu_damage_s_inv,
            mu_lysis_s_inv = mu_lysis_s_inv
        ),
        directional_flux = (
            diffusive_flux_kg_m2_s = diffusive_flux,
            radiotropic_drift_flux_kg_m2_s = radiotropic_drift_flux,
            J_B_total_kg_m2_s = J_B_total
        )
    )
end
