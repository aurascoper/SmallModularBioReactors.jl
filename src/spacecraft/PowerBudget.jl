# ==============================================================================
# src/spacecraft/PowerBudget.jl
# Spacecraft / Surface Power Budget & Strict Architectural Separation Between
# External Small Modular Fission Power System (Mars Kilopower Surface Bus) and
# Internal Small Modular Biological Reactor (Ambr 250 / Microfluidic Cassette)
# ==============================================================================

struct ExternalPowerSupplySpec
    supply_kind::Symbol
    designation::String
    max_electrical_power_W::Float64
    max_thermal_bus_power_W::Float64
    is_external_to_bioreactor::Bool

    function ExternalPowerSupplySpec(;
            supply_kind::Symbol = :mars_kilopower_fission_surface_bus,
            designation::AbstractString = "Kilopower 10 kWe Umbrella-Radiator Surface Power Unit (External Bus Only)",
            max_electrical_power_W::Real = 10000.0,
            max_thermal_bus_power_W::Real = 30000.0,
            is_external_to_bioreactor::Bool = true)
        supply_kind in (:mars_kilopower_fission_surface_bus, :spacecraft_solar_battery_bus, :benchtop_mains_bus) ||
            throw(ArgumentError("Unsupported supply_kind :$(supply_kind)"))
        (max_electrical_power_W > 0.0 && max_thermal_bus_power_W >= 0.0) ||
            throw(ArgumentError("Power capacities must be positive"))
        is_external_to_bioreactor ||
            throw(ArgumentError("Section 9 Boundary Violation: Fission/electrical power supply must remain strictly external to the biological reactor"))
        new(
            supply_kind,
            String(designation),
            Float64(max_electrical_power_W),
            Float64(max_thermal_bus_power_W),
            is_external_to_bioreactor
        )
    end
end

struct BioreactorPowerDemand
    reactor_designation::String
    pump_power_W::Float64
    heater_power_W::Float64
    optical_cytometry_sensors_W::Float64
    carousel_centrifuge_W::Float64
    simulate_fission_criticality_inside_bioreactor::Bool

    function BioreactorPowerDemand(;
            reactor_designation::AbstractString = "Small Modular Biological Reactor (Ambr 250 / BIM Microfluidic Cassette)",
            pump_power_W::Real = 18.0,
            heater_power_W::Real = 45.0,
            optical_cytometry_sensors_W::Real = 24.0,
            carousel_centrifuge_W::Real = 15.0,
            simulate_fission_criticality_inside_bioreactor::Bool = false)
        simulate_fission_criticality_inside_bioreactor && throw(ArgumentError(
            "Section 9 Hard Boundary Violation: SmallModularBioReactors.jl is a biological/biochemical bioprocess engineering layer, NOT a nuclear fission criticality simulator"
        ))
        (pump_power_W >= 0.0 && heater_power_W >= 0.0 && optical_cytometry_sensors_W >= 0.0 && carousel_centrifuge_W >= 0.0) ||
            throw(ArgumentError("Component power demands must be non-negative"))
        new(
            String(reactor_designation),
            Float64(pump_power_W),
            Float64(heater_power_W),
            Float64(optical_cytometry_sensors_W),
            Float64(carousel_centrifuge_W),
            false
        )
    end
end

"""
    evaluate_power_budget(supply::ExternalPowerSupplySpec, demand::BioreactorPowerDemand)

Evaluates electrical power margin and verifies that the external power bus (e.g. Mars Kilopower
fission surface unit or solar array) is strictly separated from the biological reactor state.
"""
function evaluate_power_budget(supply::ExternalPowerSupplySpec, demand::BioreactorPowerDemand)
    total_demand_W = demand.pump_power_W +
                     demand.heater_power_W +
                     demand.optical_cytometry_sensors_W +
                     demand.carousel_centrifuge_W
    margin_W = supply.max_electrical_power_W - total_demand_W
    margin_fraction = margin_W / supply.max_electrical_power_W

    return (
        budget_satisfied = (margin_W >= 0.0),
        total_demand_W = total_demand_W,
        available_electrical_W = supply.max_electrical_power_W,
        power_margin_W = margin_W,
        power_margin_fraction = margin_fraction,
        fission_supply_external_only = supply.is_external_to_bioreactor,
        bioreactor_is_not_fission_simulator = !demand.simulate_fission_criticality_inside_bioreactor
    )
end
