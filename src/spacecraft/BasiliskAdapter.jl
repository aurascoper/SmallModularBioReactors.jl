# ==============================================================================
# src/spacecraft/BasiliskAdapter.jl
# CU Boulder AVS Lab / LASP Basilisk Astrodynamics & Payload Simulation Adapter
# (arXiv:2605.12443 Reproducible GN&C, arXiv:2605.29308 Modular Astrobiological
# Experiments Payload, ESA BIM Carousel, scenarioMonteCarloSpice)
# ==============================================================================

struct BasiliskSpacecraftState
    scenario_name::String
    monte_carlo_run_index::Int
    mission_elapsed_time_s::Float64
    target_body::Symbol
    r_BN_N_m::Vector{Float64}
    v_BN_N_m_s::Vector{Float64}
    sigma_BN::Vector{Float64}
    omega_BN_B_rad_s::Vector{Float64}
    bim_carousel_rpm::Float64
    in_eclipse::Bool
    nutrient_valve_open::Bool

    function BasiliskSpacecraftState(;
            scenario_name::AbstractString = "scenarioMonteCarloSpice_BIM",
            monte_carlo_run_index::Int = 1,
            mission_elapsed_time_s::Real = 0.0,
            target_body::Symbol = :HST_LEO,
            r_BN_N_m::AbstractVector{<:Real} = [6.918137e6, 0.0, 0.0],
            v_BN_N_m_s::AbstractVector{<:Real} = [0.0, 7585.0, 0.0],
            sigma_BN::AbstractVector{<:Real} = [0.01, -0.02, 0.05],
            omega_BN_B_rad_s::AbstractVector{<:Real} = [0.0, 0.0, 0.0],
            bim_carousel_rpm::Real = 12.0,
            in_eclipse::Bool = false,
            nutrient_valve_open::Bool = true)
        isempty(strip(scenario_name)) && throw(ArgumentError("scenario_name cannot be empty"))
        monte_carlo_run_index >= 1 || throw(ArgumentError("monte_carlo_run_index must be >= 1"))
        mission_elapsed_time_s >= 0.0 || throw(ArgumentError("mission_elapsed_time_s must be >= 0"))
        (length(r_BN_N_m) == 3 && length(v_BN_N_m_s) == 3 && length(sigma_BN) == 3 && length(omega_BN_B_rad_s) == 3) ||
            throw(ArgumentError("Spacecraft position, velocity, MRP attitude sigma_BN, and angular rate must be 3-vectors"))
        bim_carousel_rpm >= 0.0 || throw(ArgumentError("bim_carousel_rpm must be >= 0"))

        new(
            String(scenario_name),
            monte_carlo_run_index,
            Float64(mission_elapsed_time_s),
            target_body,
            Float64.(r_BN_N_m),
            Float64.(v_BN_N_m_s),
            Float64.(sigma_BN),
            Float64.(omega_BN_B_rad_s),
            Float64(bim_carousel_rpm),
            in_eclipse,
            nutrient_valve_open
        )
    end
end

"""
    step_basilisk_spacecraft_state(state::BasiliskSpacecraftState, dt_s::Real;
        et0_s::Real = 8.46892869184e8,
        carousel_radius_m::Real = 0.14)

Advances the Basilisk spacecraft state by `dt_s` using `spice_spkpos` and computes the
ESA Biology in Microgravity (`BIM`) carousel artificial-gravity vector `g_eff` at `carousel_radius_m`.
"""
function step_basilisk_spacecraft_state(
        state::BasiliskSpacecraftState,
        dt_s::Real;
        et0_s::Real = 8.46892869184e8,
        carousel_radius_m::Real = 0.14)
    dt_s > 0.0 || throw(ArgumentError("dt_s must be > 0"))
    carousel_radius_m > 0.0 || throw(ArgumentError("carousel_radius_m must be > 0"))

    new_met_s = state.mission_elapsed_time_s + Float64(dt_s)
    et_now = Float64(et0_s) + new_met_s
    spice_out = spice_spkpos(
        state.target_body,
        et_now;
        monte_carlo_run_index = state.monte_carlo_run_index
    )

    omega_carousel_rad_s = (state.bim_carousel_rpm * 2.0 * pi) / 60.0
    total_omega = [
        state.omega_BN_B_rad_s[1],
        state.omega_BN_B_rad_s[2],
        state.omega_BN_B_rad_s[3] + omega_carousel_rad_s
    ]

    grav_body = state.target_body == :HST_LEO ? :MICROGRAVITY_LEO : state.target_body
    g_out = compute_effective_gravity(
        body = grav_body,
        omega_rad_s = total_omega,
        r_offset_m = [Float64(carousel_radius_m), 0.0, 0.0]
    )

    in_eclipse_now = (spice_out.lighting_angle_deg > 125.0)
    next_state = BasiliskSpacecraftState(
        scenario_name = state.scenario_name,
        monte_carlo_run_index = state.monte_carlo_run_index,
        mission_elapsed_time_s = new_met_s,
        target_body = state.target_body,
        r_BN_N_m = spice_out.position_si_m,
        v_BN_N_m_s = spice_out.velocity_si_m_s,
        sigma_BN = state.sigma_BN,
        omega_BN_B_rad_s = state.omega_BN_B_rad_s,
        bim_carousel_rpm = state.bim_carousel_rpm,
        in_eclipse = in_eclipse_now,
        nutrient_valve_open = state.nutrient_valve_open
    )

    return (
        spacecraft_state = next_state,
        spice_ephemeris = spice_out,
        effective_gravity = g_out,
        carousel_omega_rad_s = omega_carousel_rad_s
    )
end
