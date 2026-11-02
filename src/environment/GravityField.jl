# ==============================================================================
# src/environment/GravityField.jl
# Effective Gravity Vector g_eff with Planetary Gravity, Spacecraft Thrust/Drag
# Acceleration, Euler Angular Acceleration, Centrifugal (ESA BIM Carousel), and
# Coriolis Terms
# ==============================================================================

const CANONICAL_SURFACE_GRAVITY_M_S2 = Dict{Symbol,Float64}(
    :EARTH => 9.80665,
    :MARS => 3.72076,
    :MOON => 1.62,
    :ENCELADUS => 0.113,
    :MICROGRAVITY_LEO => 1.0e-6,
)

"""
    compute_effective_gravity(;
        body::Symbol = :MARS,
        g_planetary_m_s2::Union{Nothing,AbstractVector{<:Real}} = nothing,
        a_spacecraft_m_s2::AbstractVector{<:Real} = [0.0, 0.0, 0.0],
        omega_rad_s::AbstractVector{<:Real} = [0.0, 0.0, 0.0],
        omega_dot_rad_s2::AbstractVector{<:Real} = [0.0, 0.0, 0.0],
        r_offset_m::AbstractVector{<:Real} = [0.15, 0.0, 0.0],
        v_rel_m_s::AbstractVector{<:Real} = [0.0, 0.0, 0.0])

Computes the effective gravity vector inside a spacecraft or surface bioreactor module:
`g_eff = g_planetary - a_sc - omega_dot x r - omega x (omega x r) - 2 * omega x v_rel`
Blocks silent reuse of terrestrial `9.80665 m/s^2` when `body != :EARTH`.
"""
function compute_effective_gravity(;
        body::Symbol = :MARS,
        g_planetary_m_s2::Union{Nothing,AbstractVector{<:Real}} = nothing,
        a_spacecraft_m_s2::AbstractVector{<:Real} = [0.0, 0.0, 0.0],
        omega_rad_s::AbstractVector{<:Real} = [0.0, 0.0, 0.0],
        omega_dot_rad_s2::AbstractVector{<:Real} = [0.0, 0.0, 0.0],
        r_offset_m::AbstractVector{<:Real} = [0.15, 0.0, 0.0],
        v_rel_m_s::AbstractVector{<:Real} = [0.0, 0.0, 0.0])
    haskey(CANONICAL_SURFACE_GRAVITY_M_S2, body) ||
        throw(ArgumentError("Unsupported body :$(body); must be one of $(collect(keys(CANONICAL_SURFACE_GRAVITY_M_S2)))"))
    for (nm, v) in (
            ("a_spacecraft_m_s2", a_spacecraft_m_s2),
            ("omega_rad_s", omega_rad_s),
            ("omega_dot_rad_s2", omega_dot_rad_s2),
            ("r_offset_m", r_offset_m),
            ("v_rel_m_s", v_rel_m_s))
        length(v) == 3 || throw(ArgumentError("$(nm) must be a 3-element vector"))
        all(isfinite, v) || throw(ArgumentError("$(nm) entries must be finite"))
    end

    g_plan = if g_planetary_m_s2 === nothing
        [0.0, 0.0, -CANONICAL_SURFACE_GRAVITY_M_S2[body]]
    else
        length(g_planetary_m_s2) == 3 || throw(ArgumentError("g_planetary_m_s2 must be a 3-element vector"))
        gp = Float64.(g_planetary_m_s2)
        if body != :EARTH && isapprox(norm(gp), 9.80665; atol = 1.0e-2)
            throw(ArgumentError("Terrestrial gravity reuse blocked for non-terrestrial body :$(body)"))
        end
        gp
    end

    a_sc = Float64.(a_spacecraft_m_s2)
    w = Float64.(omega_rad_s)
    w_dot = Float64.(omega_dot_rad_s2)
    r = Float64.(r_offset_m)
    v_rel = Float64.(v_rel_m_s)

    euler_term = cross(w_dot, r)
    centrifugal_term = -cross(w, cross(w, r))
    coriolis_term = -2.0 .* cross(w, v_rel)
    g_eff = g_plan .- a_sc .- euler_term .+ centrifugal_term .+ coriolis_term

    return (
        body = body,
        g_planetary_m_s2 = g_plan,
        centrifugal_accel_m_s2 = centrifugal_term,
        coriolis_accel_m_s2 = coriolis_term,
        euler_accel_m_s2 = -euler_term,
        g_eff_m_s2 = g_eff,
        g_eff_magnitude_m_s2 = norm(g_eff),
        equivalent_earth_g = norm(g_eff) / 9.80665
    )
end
