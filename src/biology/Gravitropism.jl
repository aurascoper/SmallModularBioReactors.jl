# ==============================================================================
# src/biology/Gravitropism.jl
# Gravitropic / Sedimentation / Centrifugal Directed Biomass Flux:
#   J_B = -D_B * grad(B) + chi_g * B * g_eff
# ==============================================================================

"""
    compute_gravitropic_flux(;
        B_kg_m3::Real = 0.45,
        grad_B_kg_m4::AbstractVector{<:Real} = [0.0, 0.0, -5.0],
        g_eff_m_s2::AbstractVector{<:Real} = [0.0, 0.0, -3.72076],
        D_B_m2_s::Real = 2.0e-11,
        chi_g_s::Real = 1.2e-7)

Computes the 3D microbial/hyphal gravitropic & buoyancy-drift flux vector:
`J_B = -D_B * grad(B) + chi_g * B * g_eff` (`kg / (m^2 * s)`)
in microgravity, Martian/Lunar partial gravity, or ESA `BIM` centrifuge frames.
"""
function compute_gravitropic_flux(;
        B_kg_m3::Real = 0.45,
        grad_B_kg_m4::AbstractVector{<:Real} = [0.0, 0.0, -5.0],
        g_eff_m_s2::AbstractVector{<:Real} = [0.0, 0.0, -3.72076],
        D_B_m2_s::Real = 2.0e-11,
        chi_g_s::Real = 1.2e-7)
    B_kg_m3 >= 0.0 || throw(ArgumentError("B_kg_m3 must be >= 0"))
    D_B_m2_s > 0.0 || throw(ArgumentError("D_B_m2_s must be > 0"))
    (length(grad_B_kg_m4) == 3 && length(g_eff_m_s2) == 3) ||
        throw(ArgumentError("grad_B_kg_m4 and g_eff_m_s2 must be 3-vectors"))
    (all(isfinite, grad_B_kg_m4) && all(isfinite, g_eff_m_s2)) ||
        throw(ArgumentError("Vectors must be finite"))

    diff_term = -Float64(D_B_m2_s) .* Float64.(grad_B_kg_m4)
    grav_term = Float64(chi_g_s) * Float64(B_kg_m3) .* Float64.(g_eff_m_s2)
    J_B = diff_term .+ grav_term

    return (
        diffusive_flux_kg_m2_s = diff_term,
        gravitropic_drift_flux_kg_m2_s = grav_term,
        J_B_total_kg_m2_s = J_B,
        flux_magnitude_kg_m2_s = norm(J_B)
    )
end
