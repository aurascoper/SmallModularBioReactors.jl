# ==============================================================================
# src/spacecraft/ThermalControl.jl
# Bioreactor Payload Thermal Control & Spatial Temperature-Gradient Audit
# ==============================================================================

"""
    step_reactor_thermal_balance(;
        T_current_K::AbstractVector{<:Real} = fill(298.15, 4),
        dt_s::Real = 60.0,
        heat_capacity_J_m3_K::Real = 3.8e6,
        heater_power_W_m3::Real = 120.0,
        metabolic_heat_W_m3::Real = 15.0,
        radiolytic_heat_W_m3::Real = 0.5,
        bus_sink_temp_K::Real = 297.65,
        conductance_W_m3_K::Real = 270.0,
        cell_spacing_m::Real = 1.0e-3)

Integrates cell-scale thermal energy conservation over `dt_s` and computes the internal
spatial temperature gradient `|grad T|` (`K/m`) so thermal artifacts can be gated in
biological tropism experiments.
"""
function step_reactor_thermal_balance(;
        T_current_K::AbstractVector{<:Real} = fill(298.15, 4),
        dt_s::Real = 60.0,
        heat_capacity_J_m3_K::Real = 3.8e6,
        heater_power_W_m3::Real = 120.0,
        metabolic_heat_W_m3::Real = 15.0,
        radiolytic_heat_W_m3::Real = 0.5,
        bus_sink_temp_K::Real = 297.65,
        conductance_W_m3_K::Real = 270.0,
        cell_spacing_m::Real = 1.0e-3)
    n = length(T_current_K)
    n >= 1 || throw(ArgumentError("T_current_K must not be empty"))
    all(t -> isfinite(t) && t > 0.0, T_current_K) || throw(ArgumentError("Temperatures must be > 0 K"))
    (dt_s > 0.0 && heat_capacity_J_m3_K > 0.0 && conductance_W_m3_K >= 0.0 && cell_spacing_m > 0.0) ||
        throw(ArgumentError("dt_s, heat_capacity_J_m3_K, and cell_spacing_m must be > 0"))

    q_source = Float64(heater_power_W_m3) + Float64(metabolic_heat_W_m3) + Float64(radiolytic_heat_W_m3)
    T_next_K = zeros(Float64, n)
    for i in 1:n
        T_i = Float64(T_current_K[i])
        q_loss = Float64(conductance_W_m3_K) * (T_i - Float64(bus_sink_temp_K))
        dT_dt = (q_source - q_loss) / Float64(heat_capacity_J_m3_K)
        T_next_K[i] = T_i + dT_dt * Float64(dt_s)
    end

    max_grad_T_K_m = if n == 1
        0.0
    else
        maximum(abs(T_next_K[i + 1] - T_next_K[i]) / Float64(cell_spacing_m) for i in 1:(n - 1))
    end

    within_biocompatible_window = all(277.15 <= t <= 315.15 for t in T_next_K)
    return (
        T_next_K = T_next_K,
        mean_T_K = mean(T_next_K),
        max_grad_T_K_m = max_grad_T_K_m,
        within_biocompatible_window = within_biocompatible_window,
        thermal_gradient_negligible_for_tropism = (max_grad_T_K_m <= 5.0)
    )
end
