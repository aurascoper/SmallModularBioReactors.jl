# ==============================================================================
# src/biology/SyntrophicBioenergetics.jl
# Syntrophic Bioenergetics Workstream (SYN0..SYN6):
# - SYN0: Stoichiometric & electron-equivalent ledger (alt acceptors, >390 ppm H2 butanol shift)
# - SYN1: Two-population ODE with mandatory partner-specific qPCR biomass & Monod vs Moser PHF-HUR (70 C cutoff)
# - SYN2: Gas/liquid H2 transfer & k_L*a vs PHF-HUR K_S bias quantification
# - SYN3: Transformed Gibbs free energy Delta_r G', NutMEG closure (<280 K gate), & regulatory silence separation
# - SYN4: Fisher identifiability rank check, Bayesian calibration, & Hydrogen Hypothesis / Asgard :Incomplete quarantine
# - SYN5: Spatial aggregate Fickian diffusion (5.30 -> 0.29 um) refusing uncalibrated CPM lattice units
# - SYN6: Deterministic JSON field-bundle export for MATLAB R2026b 3D EDA
# ==============================================================================

const SYNTROPHY_MILESTONES_S0_TO_S6 = (
    :SYN0_stoichiometric_ledger,
    :SYN1_two_population_ode_and_phf_hur_kinetics,
    :SYN2_gas_liquid_transfer_and_ks_bias,
    :SYN3_thermodynamics_nutmeg_and_regulatory_silence,
    :SYN4_bayesian_identifiability_and_evolutionary_quarantine,
    :SYN5_spatial_aggregate_fickian_diffusion,
    :SYN6_matlab_3d_eda_bundle_export,
)

const QUARANTINED_SYNTROPHIC_EVOLUTIONARY_CLAIMS = (
    :hydrogen_hypothesis_ancestral_engulfment,
    :gene_transfer_nucleus_origin,
    :asgard_heimdallarchaeia_obligate_h2_syntrophy,
)

function _json_escape(s::AbstractString)
    buf = IOBuffer()
    for ch in s
        if ch == '"'
            write(buf, "\\\"")
        elseif ch == '\\'
            write(buf, "\\\\")
        elseif ch == '\n'
            write(buf, "\\n")
        elseif ch == '\r'
            write(buf, "\\r")
        elseif ch == '\t'
            write(buf, "\\t")
        else
            write(buf, ch)
        end
    end
    return String(take!(buf))
end

function _to_json_str(x; indent::Int = 0)
    pad = "  "^indent
    pad_in = "  "^(indent + 1)
    if x === nothing
        return "null"
    elseif x isa Bool
        return x ? "true" : "false"
    elseif x isa Integer
        return string(x)
    elseif x isa AbstractFloat
        isfinite(x) || throw(ArgumentError("Non-finite float cannot be serialized to strict JSON: $(x)"))
        return repr(Float64(x))
    elseif x isa Symbol
        return "\"" * _json_escape(String(x)) * "\""
    elseif x isa AbstractString
        return "\"" * _json_escape(x) * "\""
    elseif x isa AbstractVector || x isa Tuple
        isempty(x) && return "[]"
        if all(el -> el isa Number || el isa Bool, x)
            return "[" * join([_to_json_str(el; indent = 0) for el in x], ", ") * "]"
        else
            items = [pad_in * _to_json_str(el; indent = indent + 1) for el in x]
            return "[\n" * join(items, ",\n") * "\n" * pad * "]"
        end
    elseif x isa NamedTuple
        ks = keys(x)
        isempty(ks) && return "{}"
        pairs_str = [pad_in * "\"" * _json_escape(String(k)) * "\": " * _to_json_str(getfield(x, k); indent = indent + 1) for k in ks]
        return "{\n" * join(pairs_str, ",\n") * "\n" * pad * "}"
    elseif x isa AbstractDict
        ks = sort!(collect(keys(x)); by = string)
        isempty(ks) && return "{}"
        pairs_str = [pad_in * "\"" * _json_escape(string(k)) * "\": " * _to_json_str(x[k]; indent = indent + 1) for k in ks]
        return "{\n" * join(pairs_str, ",\n") * "\n" * pad * "}"
    else
        throw(ArgumentError("Unsupported JSON serialization type: $(typeof(x))"))
    end
end

"""
    evaluate_syntrophic_stoichiometric_ledger(;
        substrate::Symbol = :propionate,
        donor_consumed_mol::Real = 0.0203,
        acetate_produced_mol::Union{Nothing,Real} = nothing,
        ch4_produced_mol::Union{Nothing,Real} = nothing,
        h2_accumulated_mol::Real = 0.0,
        formate_accumulated_mol::Real = 0.0,
        alt_acceptor_so4_reduced_mol::Real = 0.0,
        alt_acceptor_no3_reduced_mol::Real = 0.0,
        p_h2_ppm::Real = 250.0,
        bicarbonate_limited::Bool = false,
        butanol_produced_mol::Real = 0.0,
        force_unbalanced_electron_ledger::Bool = false,
        rel_tol::Real = 1.0e-6)

SYN0: Verifies elemental and electron-equivalent conservation across syntrophic donor
oxidation, hydrogen/formate exchange, hydrogenotrophic methanogenesis (`CO2 + 4 H2 -> CH4 + 2 H2O`),
alternative electron acceptors (`SO4^2-`, `NO3-`), and the `>390 ppm` H2 butanol-formation shift
under bicarbonate limitation (Junicke 2015).
"""
function evaluate_syntrophic_stoichiometric_ledger(;
        substrate::Symbol = :propionate,
        donor_consumed_mol::Real = 0.0203,
        acetate_produced_mol::Union{Nothing,Real} = nothing,
        ch4_produced_mol::Union{Nothing,Real} = nothing,
        h2_accumulated_mol::Real = 0.0,
        formate_accumulated_mol::Real = 0.0,
        alt_acceptor_so4_reduced_mol::Real = 0.0,
        alt_acceptor_no3_reduced_mol::Real = 0.0,
        p_h2_ppm::Real = 250.0,
        bicarbonate_limited::Bool = false,
        butanol_produced_mol::Real = 0.0,
        force_unbalanced_electron_ledger::Bool = false,
        rel_tol::Real = 1.0e-6)
    force_unbalanced_electron_ledger && throw(ArgumentError(
        "SYN0 Stoichiometric Violation: unbalanced electron-donor ledger refused"
    ))
    substrate in (:propionate, :lactate, :butyrate, :glucose) ||
        throw(ArgumentError("Unsupported substrate :$(substrate); expected :propionate, :lactate, :butyrate, or :glucose"))
    donor_m = Float64(donor_consumed_mol)
    h2_m = Float64(h2_accumulated_mol)
    fmt_m = Float64(formate_accumulated_mol)
    so4_m = Float64(alt_acceptor_so4_reduced_mol)
    no3_m = Float64(alt_acceptor_no3_reduced_mol)
    but_m = Float64(butanol_produced_mol)
    p_h2 = Float64(p_h2_ppm)
    tol = Float64(rel_tol)

    (donor_m >= 0.0 && h2_m >= 0.0 && fmt_m >= 0.0 && so4_m >= 0.0 && no3_m >= 0.0 && but_m >= 0.0 && p_h2 >= 0.0 && tol > 0.0) ||
        throw(ArgumentError("Stoichiometric molar amounts, p_h2_ppm, and rel_tol must be non-negative"))

    # Stoichiometric yields per mole donor oxidized to acetate + HCO3- + H2
    # Each H2 (or formate) carries 2 e- equivalents.
    h2_equiv_per_donor, ac_per_donor, ic_per_donor = if substrate == :propionate
        (3.0, 1.0, 1.0)
    elseif substrate == :lactate
        (2.0, 1.0, 1.0)
    elseif substrate == :butyrate
        (2.0, 2.0, 0.0)
    else # :glucose
        (4.0, 2.0, 2.0)
    end

    e_donor_released_mol = 2.0 * h2_equiv_per_donor * donor_m

    # Sinks in electron-moles:
    # CH4 from CO2 + 4 H2 -> CH4 + 2 H2O consumes 8 e- per mol CH4
    # H2 consumes 2 e- per mol; Formate consumes 2 e- per mol
    # SO4^2- -> H2S consumes 8 e- per mol; NO3- -> NO2- consumes 2 e- per mol
    # 2 Acetate + 4 H2 + 2 H+ -> Butanol + 3 H2O consumes 8 e- per mol butanol (and 2 mol acetate)
    e_non_ch4 = 2.0 * h2_m + 2.0 * fmt_m + 8.0 * so4_m + 2.0 * no3_m + 8.0 * but_m
    ch4_m = ch4_produced_mol === nothing ? max(0.0, (e_donor_released_mol - e_non_ch4) / 8.0) : Float64(ch4_produced_mol)
    ch4_m >= 0.0 || throw(ArgumentError("ch4_produced_mol cannot be negative"))

    ac_expected_mol = ac_per_donor * donor_m - 2.0 * but_m
    if acetate_produced_mol === nothing
        ac_expected_mol >= -1.0e-15 || throw(ArgumentError(
            "SYN0 Stoichiometric Violation: butanol_produced_mol ($(but_m)) exceeds available acetyl-donor stoichiometry without external acetate input"
        ))
        ac_m = max(0.0, ac_expected_mol)
    else
        ac_m = Float64(acetate_produced_mol)
    end
    ac_m >= 0.0 || throw(ArgumentError("acetate_produced_mol cannot be negative"))

    e_ch4_mol = 8.0 * ch4_m
    e_total_accounted_mol = e_ch4_mol + e_non_ch4
    e_abs_residual = abs(e_donor_released_mol - e_total_accounted_mol)
    e_rel_residual = e_abs_residual / max(1.0e-18, max(e_donor_released_mol, e_total_accounted_mol))

    # Independent carbon balance:
    # Donor oxidation releases (2*ac_per_donor + ic_per_donor)*donor_m C.
    # Methanogenesis (CO2 + 4 H2 -> CH4) and formate accumulation each draw 1 mol C from inorganic C.
    # 2 Acetate -> Butanol preserves reduced organic carbon.
    c_per_donor = substrate == :propionate ? 3.0 : (substrate == :lactate ? 3.0 : (substrate == :butyrate ? 4.0 : 6.0))
    c_donor = c_per_donor * donor_m
    c_reduced_products = 2.0 * ac_m + 4.0 * but_m + ch4_m + fmt_m
    net_inorganic_carbon_change = ic_per_donor * donor_m - ch4_m - fmt_m
    c_total_accounted = c_reduced_products + net_inorganic_carbon_change
    c_balance_rel_residual = abs(c_total_accounted - c_donor) / max(1.0e-18, max(c_donor, c_total_accounted))
    acetate_balance_rel_err = abs(ac_m - max(0.0, ac_expected_mol)) / max(1.0e-18, max(ac_m, abs(ac_expected_mol)))
    c_expected_from_donor_only = ac_per_donor * donor_m * 2.0 + ic_per_donor * donor_m

    butanol_shift_triggered = (p_h2 > 390.0) && bicarbonate_limited
    butanol_shift_active = butanol_shift_triggered || (but_m > 0.0)

    return (
        stage_id = :SYN0_stoichiometric_ledger,
        substrate = substrate,
        donor_consumed_mol = donor_m,
        acetate_produced_mol = ac_m,
        acetate_expected_mol = max(0.0, ac_expected_mol),
        acetate_balance_closed = (acetate_balance_rel_err <= tol),
        ch4_produced_mol = ch4_m,
        h2_accumulated_mol = h2_m,
        formate_accumulated_mol = fmt_m,
        alt_acceptor_so4_reduced_mol = so4_m,
        alt_acceptor_no3_reduced_mol = no3_m,
        butanol_produced_mol = but_m,
        p_h2_ppm = p_h2,
        bicarbonate_limited = bicarbonate_limited,
        butanol_shift_triggered_above_390ppm = butanol_shift_triggered,
        butanol_shift_active = butanol_shift_active,
        e_donor_released_mol = e_donor_released_mol,
        e_ch4_sink_mol = e_ch4_mol,
        e_alternative_acceptors_mol = 8.0 * so4_m + 2.0 * no3_m,
        e_total_accounted_mol = e_total_accounted_mol,
        electron_balance_abs_residual_mol = e_abs_residual,
        electron_balance_rel_residual = e_rel_residual,
        electron_balance_closed = (e_rel_residual <= tol) && (c_balance_rel_residual <= tol),
        carbon_in_mol = c_donor,
        carbon_reduced_products_mol = c_reduced_products,
        net_inorganic_carbon_mol = net_inorganic_carbon_change,
        carbon_balance_rel_residual = c_balance_rel_residual,
        carbon_balance_closed = (c_balance_rel_residual <= tol),
        reference_c_expected_mol = c_expected_from_donor_only
    )
end

function _phf_hur_kinetic_params(temperature_C::Float64)
    # Fachal-Suarez et al. (2025), Chem. Eng. J. Adv. 10.1016/j.ceja.2025.100951:
    # - K_S(27 C) = 4.27e-3 mM (= 4.27e-3 mol/m^3), K_S(65 C) = 6.87e-2 mM (= 6.87e-2 mol/m^3)
    # - Uptake rate peaks at 773 mmol H2 / (g VSS * d) at 65 C (= 8.946759e-3 mol/(kg*s))
    # - Nil (0.0) at >= 70 C (thermal inactivation cutoff)
    cutoff_70C = (temperature_C >= 70.0)
    frac_temp = clamp((temperature_C - 27.0) / (65.0 - 27.0), 0.0, 1.0)
    k_s_mM = 4.27e-3 * ((6.87e-2 / 4.27e-3)^frac_temp)
    q_65C_mol_kg_s = 773.0 / 86400.0
    if cutoff_70C
        q_max_mol_kg_s = 0.0
    elseif temperature_C <= 65.0
        # Smooth Q10-like ramp from ~160 at 27 C up to 773 mmol/(gVSS*d) at 65 C
        q_max_mol_kg_s = q_65C_mol_kg_s * (0.207^((65.0 - temperature_C) / 38.0))
    else
        # Sharp thermal collapse between 65 C and 70 C reaching 0 at 70 C
        q_max_mol_kg_s = q_65C_mol_kg_s * max(0.0, 1.0 - (temperature_C - 65.0) / 5.0)
    end
    return (
        q_max_H2_mol_kg_s = q_max_mol_kg_s,
        K_S_H2_mol_m3 = k_s_mM,
        thermal_cutoff_at_70C = cutoff_70C
    )
end

"""
    step_syntrophic_two_population_kinetics(;
        X_donor_kg_m3::Real = 0.35,
        X_methanogen_kg_m3::Real = 0.15,
        use_total_biomass_only::Bool = false,
        S_donor_mol_m3::Real = 10.0,
        H2_aq_mol_m3::Real = 2.0e-2,
        temperature_C::Real = 37.0,
        dt_s::Real = 60.0,
        n_steps::Int = 20,
        kinetic_model::Symbol = :moser,
        n_moser::Real = 1.65,
        K_I_H2_mol_m3::Real = 0.05,
        q_max_donor_mol_kg_s::Real = 8.0e-4,
        Y_donor_kg_mol::Real = 0.012,
        Y_methanogen_kg_mol_H2::Real = 0.0025,
        V_reactor_m3::Real = 1.2e-4)

SYN1: Integrates the coupled donor-methanogen ODE enforcing partner-specific qPCR biomass
(`X_donor`, `X_methanogen` separate; refuses `use_total_biomass_only = true` per Junicke 2015),
comparing Monod vs Moser kinetics across the Fachal-Suarez (2025) `PHF-HUR` temperature window
with strict zero uptake at `T >= 70 C`, and verifying electron-equivalent conservation.
"""
function step_syntrophic_two_population_kinetics(;
        X_donor_kg_m3::Real = 0.35,
        X_methanogen_kg_m3::Real = 0.15,
        use_total_biomass_only::Bool = false,
        S_donor_mol_m3::Real = 10.0,
        H2_aq_mol_m3::Real = 2.0e-2,
        temperature_C::Real = 37.0,
        dt_s::Real = 60.0,
        n_steps::Int = 20,
        kinetic_model::Symbol = :moser,
        n_moser::Real = 1.65,
        K_I_H2_mol_m3::Real = 0.05,
        q_max_donor_mol_kg_s::Real = 8.0e-4,
        Y_donor_kg_mol::Real = 0.012,
        Y_methanogen_kg_mol_H2::Real = 0.0025,
        V_reactor_m3::Real = 1.2e-4)
    use_total_biomass_only && throw(ArgumentError(
        "SYN1 Partner Biomass Violation (Junicke 2015): total biomass alone cannot identify partner-specific electron fluxes; separate qPCR X_donor and X_methanogen are required"
    ))
    kinetic_model in (:monod, :moser) ||
        throw(ArgumentError("kinetic_model must be :monod or :moser, got :$(kinetic_model)"))

    Xd = Float64(X_donor_kg_m3)
    Xm = Float64(X_methanogen_kg_m3)
    Sd = Float64(S_donor_mol_m3)
    H2 = Float64(H2_aq_mol_m3)
    Tc = Float64(temperature_C)
    dt = Float64(dt_s)
    n_mos = Float64(n_moser)
    K_I = Float64(K_I_H2_mol_m3)
    q_d_max = Float64(q_max_donor_mol_kg_s)
    Yd = Float64(Y_donor_kg_mol)
    Ym = Float64(Y_methanogen_kg_mol_H2)
    V = Float64(V_reactor_m3)

    (isfinite(Xd) && Xd >= 0.0 && isfinite(Xm) && Xm >= 0.0 && isfinite(Sd) && Sd >= 0.0 &&
     isfinite(H2) && H2 >= 0.0 && isfinite(Tc) && Tc > -50.0 && dt > 0.0 && n_steps >= 1 &&
     n_mos > 0.0 && K_I > 0.0 && q_d_max >= 0.0 && Yd >= 0.0 && Ym >= 0.0 && V > 0.0) ||
        throw(ArgumentError("Invalid non-finite or negative parameters in step_syntrophic_two_population_kinetics"))

    phf = _phf_hur_kinetic_params(Tc)
    q_m_max = phf.q_max_H2_mol_kg_s
    K_S = phf.K_S_H2_mol_m3

    # Instantaneous Monod vs Moser rate comparison at initial dissolved H2
    monod_q0 = phf.thermal_cutoff_at_70C ? 0.0 : q_m_max * (H2 / (K_S + H2))
    moser_q0 = phf.thermal_cutoff_at_70C ? 0.0 : q_m_max * ((H2^n_mos) / (K_S^n_mos + H2^n_mos))

    # Compare Monod vs Moser fits against synthetic PHF-HUR ground-truth sigmoid curve (Fachal-Suarez 2025)
    h2_probe_grid = [0.002, 0.005, 0.015, 0.040, 0.080, 0.200]
    true_n = 1.65
    y_true = [phf.thermal_cutoff_at_70C ? 0.0 : q_m_max * ((h^true_n) / (K_S^true_n + h^true_n)) for h in h2_probe_grid]
    y_monod = [phf.thermal_cutoff_at_70C ? 0.0 : q_m_max * (h / (K_S + h)) for h in h2_probe_grid]
    y_moser = [phf.thermal_cutoff_at_70C ? 0.0 : q_m_max * ((h^n_mos) / (K_S^n_mos + h^n_mos)) for h in h2_probe_grid]
    sse_monod = sum((y_monod[i] - y_true[i])^2 for i in eachindex(h2_probe_grid))
    sse_moser = sum((y_moser[i] - y_true[i])^2 for i in eachindex(h2_probe_grid))

    # Conservative sub-stepped integration of:
    #   propionate + 3 H2O -> acetate + HCO3- + H+ + 3 H2  (6 e- / mol donor)
    #   CO2 + 4 H2 -> CH4 + 2 H2O                           (8 e- / mol CH4)
    Sd_init = Sd
    H2_init = H2
    CH4 = 0.0
    substeps = 20
    dt_sub = dt / substeps
    K_Sd = 0.25

    q_e_donor_last = 0.0
    q_e_methanogen_last = 0.0

    for _step in 1:n_steps
        for _sub in 1:substeps
            r_donor = (Sd > 0.0 && !phf.thermal_cutoff_at_70C) ?
                      (q_d_max * Xd * (Sd / (K_Sd + Sd)) * (1.0 / (1.0 + H2 / K_I))) : 0.0
            dSd_max = Sd / dt_sub
            r_donor = min(r_donor, 0.90 * dSd_max)

            h2_prod_rate = 3.0 * r_donor
            q_h2_spec = if phf.thermal_cutoff_at_70C || H2 <= 0.0
                0.0
            elseif kinetic_model == :monod
                q_m_max * (H2 / (K_S + H2))
            else
                q_m_max * ((H2^n_mos) / (K_S^n_mos + H2^n_mos))
            end
            r_h2_cons = min(q_h2_spec * Xm, 0.90 * (H2 / dt_sub + h2_prod_rate))

            Sd = max(0.0, Sd - r_donor * dt_sub)
            H2 = max(0.0, H2 + (h2_prod_rate - r_h2_cons) * dt_sub)
            CH4 += 0.25 * r_h2_cons * dt_sub

            Xd += Yd * r_donor * dt_sub
            Xm += Ym * r_h2_cons * dt_sub

            q_e_donor_last = Xd > 0.0 ? (6.0 * r_donor / Xd) : 0.0
            q_e_methanogen_last = Xm > 0.0 ? (2.0 * r_h2_cons / Xm) : 0.0
        end
    end

    # Electron-equivalent balance over the reactor volume V
    # Initial e- in donor + initial H2 == Final e- in donor + final H2 + CH4 formed
    e_initial_mol = V * (6.0 * Sd_init + 2.0 * H2_init)
    e_final_mol = V * (6.0 * Sd + 2.0 * H2 + 8.0 * CH4)
    e_abs_res = abs(e_final_mol - e_initial_mol)
    e_rel_res = e_abs_res / max(1.0e-18, max(e_initial_mol, e_final_mol))

    return (
        stage_id = :SYN1_two_population_ode_and_phf_hur_kinetics,
        partner_qpcr_enforced = true,
        kinetic_model = kinetic_model,
        temperature_C = Tc,
        thermal_cutoff_at_70C = phf.thermal_cutoff_at_70C,
        q_max_H2_mol_kg_s = q_m_max,
        q_max_H2_mmol_gVSS_d = q_m_max * 86400.0,
        K_S_H2_mol_m3 = K_S,
        monod_specific_rate_mol_kg_s = monod_q0,
        moser_specific_rate_mol_kg_s = moser_q0,
        sse_monod = sse_monod,
        sse_moser = sse_moser,
        moser_preferred_over_monod = (!phf.thermal_cutoff_at_70C) && (sse_moser <= sse_monod),
        X_donor_final_kg_m3 = Xd,
        X_methanogen_final_kg_m3 = Xm,
        S_donor_final_mol_m3 = Sd,
        H2_aq_final_mol_m3 = H2,
        CH4_produced_mol_m3 = CH4,
        q_e_donor_mol_kg_s = q_e_donor_last,
        q_e_methanogen_mol_kg_s = q_e_methanogen_last,
        electron_balance_rel_residual = e_rel_res,
        electron_balance_conserved = (e_rel_res < 1.0e-12)
    )
end

"""
    evaluate_syntrophic_gas_liquid_transfer(;
        H2_aq_mol_m3::Real = 0.05,
        p_H2_atm::Real = 2.5e-4,
        temperature_K::Real = 310.15,
        k_L_a_s_inv::Real = 2.0e-3,
        q_H2_max_mol_kg_s::Real = 4.0e-3,
        X_methanogen_kg_m3::Real = 0.25,
        K_S_true_mol_m3::Real = 1.5e-2,
        headspace_free_phf_hur::Bool = false)

SYN2: Evaluates gas-liquid H2 mass transfer `J_gl = k_L*a * (H_H2*p_H2 - H2_aq)` and
quantifies the systematic `K_S,H2` bias when using off-gas measurements under mass-transfer
limitation vs headspace-free `PHF-HUR` dissolved-H2 probing (Junicke 2015; Fachal-Suarez 2025).
"""
function evaluate_syntrophic_gas_liquid_transfer(;
        H2_aq_mol_m3::Real = 0.05,
        p_H2_atm::Real = 2.5e-4,
        temperature_K::Real = 310.15,
        k_L_a_s_inv::Real = 2.0e-3,
        q_H2_max_mol_kg_s::Real = 4.0e-3,
        X_methanogen_kg_m3::Real = 0.25,
        K_S_true_mol_m3::Real = 1.5e-2,
        headspace_free_phf_hur::Bool = false)
    (H2_aq_mol_m3 >= 0.0 && p_H2_atm >= 0.0 && temperature_K > 200.0 &&
     k_L_a_s_inv > 0.0 && q_H2_max_mol_kg_s >= 0.0 && X_methanogen_kg_m3 >= 0.0 && K_S_true_mol_m3 > 0.0) ||
        throw(ArgumentError("Invalid physical parameters in evaluate_syntrophic_gas_liquid_transfer"))

    T = Float64(temperature_K)
    # Henry's law constant for H2 in water (0.78 mol/(m^3*atm) at 298.15 K, van 't Hoff enthalpy -4180 J/mol)
    R_gas = 8.314462618
    H_H2_mol_m3_atm = 0.78 * exp((-4180.0 / R_gas) * (1.0 / T - 1.0 / 298.15))
    H2_eq_mol_m3 = H_H2_mol_m3_atm * Float64(p_H2_atm)
    H2_aq = Float64(H2_aq_mol_m3)
    kla = Float64(k_L_a_s_inv)

    flux_gas_to_liq_mol_m3_s = headspace_free_phf_hur ? 0.0 : kla * (H2_eq_mol_m3 - H2_aq)
    oversaturation_ratio = H2_eq_mol_m3 > 0.0 ? (H2_aq / H2_eq_mol_m3) : Inf
    liquid_oversaturated = (H2_aq > H2_eq_mol_m3)

    # Damkohler-like ratio Da = (q_max * X) / (k_L*a * K_S_true)
    # When supplying H2 from gas under consumption, apparent K_S is overestimated by (1 + Da).
    # Under liquid-phase H2 production oversaturation, off-gas underestimates apparent K_S.
    da_transfer = headspace_free_phf_hur ? 0.0 : (Float64(q_H2_max_mol_kg_s) * Float64(X_methanogen_kg_m3)) / (kla * Float64(K_S_true_mol_m3))
    ks_apparent_gas_supplied = Float64(K_S_true_mol_m3) * (1.0 + da_transfer)
    ks_bias_factor = ks_apparent_gas_supplied / Float64(K_S_true_mol_m3)
    ks_apparent_liquid_oversat = (!headspace_free_phf_hur && liquid_oversaturated && oversaturation_ratio > 1.0) ?
                                 Float64(K_S_true_mol_m3) / oversaturation_ratio : Float64(K_S_true_mol_m3)

    return (
        stage_id = :SYN2_gas_liquid_transfer_and_ks_bias,
        headspace_free_phf_hur = headspace_free_phf_hur,
        gas_liquid_transfer_limited = !headspace_free_phf_hur && (da_transfer > 0.05 || abs(H2_aq - H2_eq_mol_m3) > 1.0e-6),
        henry_constant_H2_mol_m3_atm = H_H2_mol_m3_atm,
        H2_eq_mol_m3 = H2_eq_mol_m3,
        H2_aq_mol_m3 = H2_aq,
        k_L_a_s_inv = kla,
        flux_gas_to_liquid_mol_m3_s = flux_gas_to_liq_mol_m3_s,
        liquid_oversaturated = liquid_oversaturated,
        oversaturation_ratio = oversaturation_ratio,
        K_S_true_mol_m3 = Float64(K_S_true_mol_m3),
        K_S_apparent_gas_supplied_mol_m3 = ks_apparent_gas_supplied,
        K_S_apparent_liquid_oversat_mol_m3 = ks_apparent_liquid_oversat,
        ks_bias_ratio = ks_bias_factor,
        phf_hur_eliminates_kla_bias = headspace_free_phf_hur && (abs(ks_bias_factor - 1.0) < 1.0e-12)
    )
end

"""
    evaluate_syntrophic_thermodynamics(;
        substrate::Symbol = :propionate,
        temperature_K::Real = 310.15,
        donor_mol_L::Real = 0.010,
        acetate_mol_L::Real = 0.005,
        hco3_mol_L::Real = 0.030,
        p_H2_atm::Real = 1.0e-4,
        p_CH4_atm::Real = 0.30,
        pH::Real = 7.0,
        equate_thermodynamic_feasibility_with_monoculture_flux::Bool = false)

SYN3: Computes transformed reaction Gibbs free energies `Delta_r G'` for syntrophic
donor oxidation and hydrogenotrophic methanogenesis, connects to `NutMEG` maintenance
budget closure (`<280 K` quarantine), and enforces regulatory silence separation
(Hamilton et al. 2015: thermodynamic feasibility `Delta_r G' < 0` alone does not
explain absence of H2 production in *S. fumaroxidans* monoculture).
"""
function evaluate_syntrophic_thermodynamics(;
        substrate::Symbol = :propionate,
        temperature_K::Real = 310.15,
        donor_mol_L::Real = 0.010,
        acetate_mol_L::Real = 0.005,
        hco3_mol_L::Real = 0.030,
        p_H2_atm::Real = 1.0e-4,
        p_CH4_atm::Real = 0.30,
        pH::Real = 7.0,
        equate_thermodynamic_feasibility_with_monoculture_flux::Bool = false)
    equate_thermodynamic_feasibility_with_monoculture_flux && throw(ArgumentError(
        "SYN3 Regulatory Silence Violation (Hamilton et al. 2015): thermodynamic feasibility (Delta_r G' < 0) cannot be equated with monoculture H2 flux without cellular regulatory constraints"
    ))
    substrate in (:propionate, :butyrate, :lactate) ||
        throw(ArgumentError("Unsupported substrate :$(substrate) in evaluate_syntrophic_thermodynamics"))
    (temperature_K > 0.0 && donor_mol_L > 0.0 && acetate_mol_L > 0.0 && hco3_mol_L > 0.0 &&
     p_H2_atm > 0.0 && p_CH4_atm > 0.0 && 0.0 <= pH <= 14.0) ||
        throw(ArgumentError("Thermodynamic concentrations, partial pressures, and T must be positive"))

    T = Float64(temperature_K)
    R_gas = 8.314462618
    # pH correction relative to pH 7 reference
    h_rel_pH7 = 10.0^(-(Float64(pH) - 7.0))

    delta_G0_prime_donor_J_mol, q_donor = if substrate == :propionate
        # Propionate- + 3 H2O -> Acetate- + HCO3- + H+ + 3 H2 (Delta G0' = +76,100 J/mol)
        (76100.0, (Float64(acetate_mol_L) * Float64(hco3_mol_L) * h_rel_pH7 * (Float64(p_H2_atm)^3)) / Float64(donor_mol_L))
    elseif substrate == :butyrate
        # Butyrate- + 2 H2O -> 2 Acetate- + H+ + 2 H2 (Delta G0' = +48,100 J/mol)
        (48100.0, ((Float64(acetate_mol_L)^2) * h_rel_pH7 * (Float64(p_H2_atm)^2)) / Float64(donor_mol_L))
    else
        # Lactate- + 2 H2O -> Acetate- + HCO3- + H+ + 2 H2 (Delta G0' = -4,200 J/mol)
        (-4200.0, (Float64(acetate_mol_L) * Float64(hco3_mol_L) * h_rel_pH7 * (Float64(p_H2_atm)^2)) / Float64(donor_mol_L))
    end
    delta_G_donor_J_mol = delta_G0_prime_donor_J_mol + R_gas * T * log(max(1.0e-30, q_donor))

    # Hydrogenotrophic methanogenesis: HCO3- + H+ + 4 H2 -> CH4 + 3 H2O (Delta G0' = -135,600 J/mol)
    q_meth = Float64(p_CH4_atm) / (Float64(hco3_mol_L) * h_rel_pH7 * (Float64(p_H2_atm)^4))
    delta_G_methanogen_J_mol = -135600.0 + R_gas * T * log(max(1.0e-30, q_meth))

    # Thermodynamic potential factor F_T (Jin & Bethke 2003) with minimum quantum Delta G_ATP ~ 20 kJ/mol
    delta_G_atp_quantum = 20000.0
    F_T_donor = delta_G_donor_J_mol < -delta_G_atp_quantum ?
                max(0.0, 1.0 - exp((delta_G_donor_J_mol + delta_G_atp_quantum) / (R_gas * T))) : 0.0
    F_T_meth = delta_G_methanogen_J_mol < -delta_G_atp_quantum ?
               max(0.0, 1.0 - exp((delta_G_methanogen_J_mol + delta_G_atp_quantum) / (2.0 * R_gas * T))) : 0.0

    # Connect to NutMEG maintenance closure
    dissolved_h2_mol_L = (0.78e-3) * Float64(p_H2_atm)
    nutmeg = evaluate_nutmeg_methanogen_budget_closure(
        temperature_K = T,
        dissolved_h2_mol_L = dissolved_h2_mol_L,
        target_world = :EARTH
    )

    donor_exergonic = (delta_G_donor_J_mol < 0.0)
    meth_exergonic = (delta_G_methanogen_J_mol < 0.0)
    both_exergonic = donor_exergonic && meth_exergonic
    verdict = nutmeg.out_of_calibrated_maintenance_range_below_280K ? :Incomplete : (both_exergonic ? :Accepted : :Incomplete)

    return (
        stage_id = :SYN3_thermodynamics_nutmeg_and_regulatory_silence,
        substrate = substrate,
        temperature_K = T,
        p_H2_atm = Float64(p_H2_atm),
        delta_r_G_donor_J_mol = delta_G_donor_J_mol,
        delta_r_G_methanogen_J_mol = delta_G_methanogen_J_mol,
        thermodynamic_factor_F_T_donor = F_T_donor,
        thermodynamic_factor_F_T_methanogen = F_T_meth,
        donor_oxidation_exergonic = donor_exergonic,
        methanogenesis_exergonic = meth_exergonic,
        syntrophy_thermodynamically_feasible = both_exergonic,
        thermodynamic_feasibility_sufficient_for_monoculture_h2 = false,
        regulatory_constraints_required_note = "Hamilton et al. 2015: S. fumaroxidans monoculture suppresses H2 even when exergonic; regulation/directionality required",
        nutmeg_closure = nutmeg,
        verdict = verdict
    )
end

"""
    evaluate_syntrophic_identifiability_and_calibration(;
        has_qpcr_partner_biomass::Bool = true,
        has_dissolved_h2_probe::Bool = true,
        headspace_free_phf_hur::Bool = true,
        n_mcmc_steps::Int = 160,
        seed::Int = 42,
        promote_evolutionary_narrative_to_accepted::Bool = false)

SYN4: Evaluates Fisher Information Matrix rank/conditioning for kinetic/transport
identifiability (`[q_max_donor, q_max_methanogen, K_S_H2, n_moser, k_L_a]`), runs
Hastings-Metropolis MCMC calibration on `PHF-HUR` dissolved H2 kinetics, and
quarantines Hydrogen Hypothesis / Asgard evolutionary narratives in `:Incomplete`.
"""
function evaluate_syntrophic_identifiability_and_calibration(;
        has_qpcr_partner_biomass::Bool = true,
        has_dissolved_h2_probe::Bool = true,
        headspace_free_phf_hur::Bool = true,
        n_mcmc_steps::Int = 160,
        seed::Int = 42,
        promote_evolutionary_narrative_to_accepted::Bool = false)
    promote_evolutionary_narrative_to_accepted && throw(ArgumentError(
        "SYN4 Epistemic Boundary Violation: Hydrogen Hypothesis engulfment, nucleus origin, and Asgard/Heimdallarchaeia enrichment narratives remain :Incomplete and cannot be promoted to :Accepted"
    ))

    # Build sensitivity Jacobian J (rows: observable timepoints/channels, cols: parameters)
    # Parameters:
    #   If headspace_free_phf_hur == true: theta = [q_max_donor, q_max_methanogen, K_S_H2, n_moser] (4 params)
    #   If headspace_free_phf_hur == false: theta = [q_max_donor, q_max_methanogen, K_S_H2, n_moser, k_L_a] (5 params)
    h2_levels = [0.004, 0.015, 0.045, 0.120]
    rows = Vector{Vector{Float64}}()
    for h in h2_levels
        # Channel 1: dissolved H2 uptake curve (or gas-phase CH4 rate)
        # Without dissolved H2 probe, K_S_H2 and k_L_a become collinear and Moser n_moser collapses
        if has_dissolved_h2_probe
            d_qd = has_qpcr_partner_biomass ? 0.35 * (1.0 / (1.0 + h / 0.05)) : 0.50 * (1.0 / (1.0 + h / 0.05))
            d_qm = has_qpcr_partner_biomass ? 0.15 * ((h^1.65) / (0.02^1.65 + h^1.65)) : 0.50 * (1.0 / (1.0 + h / 0.05)) # collinear without qPCR!
            d_ks = -0.15 * 1.65 * (h^1.65) * (0.02^0.65) / ((0.02^1.65 + h^1.65)^2)
            d_n = 0.15 * ((h^1.65) * (0.02^1.65) * log(max(1.0e-9, h / 0.02))) / ((0.02^1.65 + h^1.65)^2)
            if headspace_free_phf_hur
                push!(rows, [d_qd, d_qm, d_ks, d_n])
            else
                d_kla = 0.4 * d_ks + 0.02 * (0.05 - h)
                push!(rows, [d_qd, d_qm, d_ks, d_n, d_kla])
            end
        else
            # Only bulk gas CH4 observed: K_S and k_L_a collinear, n_moser unidentifiable
            d_bulk = 0.50 / (1.0 + h / 0.05)
            if headspace_free_phf_hur
                push!(rows, [d_bulk, d_bulk, 0.1 * d_bulk, 0.0])
            else
                push!(rows, [d_bulk, d_bulk, 0.1 * d_bulk, 0.0, 0.2 * d_bulk])
            end
        end
    end
    if has_qpcr_partner_biomass
        # Direct qPCR partner biomass observations separate q_max_donor from q_max_methanogen
        if headspace_free_phf_hur
            push!(rows, [1.0, 0.0, 0.0, 0.0])
            push!(rows, [0.0, 1.0, 0.0, 0.0])
        else
            push!(rows, [1.0, 0.0, 0.0, 0.0, 0.0])
            push!(rows, [0.0, 1.0, 0.0, 0.0, 0.0])
            push!(rows, [0.0, 0.0, 0.0, 0.0, 1.0]) # independent gas-transfer tracer
        end
    end

    J = reduce(vcat, [r' for r in rows])
    F = J' * J
    sv = svdvals(F)
    tol_rank = 1.0e-8 * max(1.0, maximum(sv))
    fisher_rank = count(s -> s > tol_rank, sv)
    n_params = size(J, 2)
    full_rank = (fisher_rank == n_params)
    cond_num = minimum(sv) > 1.0e-15 ? (maximum(sv) / minimum(sv)) : Inf

    # Bayesian calibration via HastingsMetropolis over [q_max_H2_mmol_gVSS_d, K_S_mM, n_moser]
    prior = GaussianPriorSpec(
        param_names = [:q_max_H2_mmol_gVSS_d, :K_S_H2_mM, :n_moser],
        prior_mean = [750.0, 0.065, 1.60],
        prior_std = [60.0, 0.012, 0.20]
    )
    h2_cal = [0.01, 0.05, 0.15]
    forward_phf = function (theta::AbstractVector{<:Real})
        qm = max(10.0, Float64(theta[1]))
        ks = max(1.0e-4, Float64(theta[2]))
        nm = clamp(Float64(theta[3]), 0.5, 4.0)
        return [qm * ((h^nm) / (ks^nm + h^nm)) / 773.0 for h in h2_cal]
    end
    y_obs = forward_phf([773.0, 0.0687, 1.65])
    mcmc = run_hastings_metropolis_calibration(
        y_obs = y_obs,
        prior = prior,
        forward_pred_fn = forward_phf,
        obs_noise_std = 0.03,
        proposal_scales = [18.0, 0.004, 0.06],
        n_steps = n_mcmc_steps,
        seed = seed
    )

    evolutionary_quarantine = (
        hydrogen_hypothesis_ancestral_engulfment = evaluate_neuroarc_evidence_gate(
            claim_id = "Martin & Muller 1998 Hydrogen Hypothesis ancestral archaeal-bacterial engulfment sequence",
            has_witness_evidence = false
        ),
        gene_transfer_nucleus_origin = evaluate_neuroarc_evidence_gate(
            claim_id = "Endosymbiotic gene transfer and nucleus-origin sequence inferred solely from modern syntrophy",
            has_witness_evidence = false
        ),
        asgard_heimdallarchaeia_obligate_h2_syntrophy = evaluate_neuroarc_evidence_gate(
            claim_id = "Asgard / Heimdallarchaeia enrichment proves exact ancestral H2-dependent partnership",
            has_witness_evidence = false
        ),
    )

    identifiability_verdict = full_rank ? :Accepted : :Incomplete

    return (
        stage_id = :SYN4_bayesian_identifiability_and_evolutionary_quarantine,
        has_qpcr_partner_biomass = has_qpcr_partner_biomass,
        has_dissolved_h2_probe = has_dissolved_h2_probe,
        headspace_free_phf_hur = headspace_free_phf_hur,
        n_parameters = n_params,
        fisher_rank = fisher_rank,
        identifiable = full_rank,
        fisher_condition_number = cond_num,
        identifiability_verdict = identifiability_verdict,
        mcmc_posterior = mcmc,
        evolutionary_quarantine = evolutionary_quarantine,
        evolutionary_track_separated = all(g.verdict == :Incomplete && !g.eligible_for_memory_banking for g in values(evolutionary_quarantine))
    )
end

"""
    evaluate_syntrophic_spatial_diffusion(;
        interspecies_distance_um::Real = 0.29,
        coordinate_unit::Symbol = :micrometers,
        lattice_pitch_um::Union{Nothing,Real} = nothing,
        D_H2_m2_s::Real = 4.5e-9,
        delta_C_H2_mol_m3::Real = 5.0e-4,
        aggregate_volume_fraction::Real = 0.15,
        specific_exchange_area_m_inv::Real = 1.5e5)

SYN5: Evaluates Fickian interspecies H2 diffusive flux across syntrophic coaggregates
(Felchner-Zwirello et al. 2013: interbacterial distance shrinking `5.30 -> 0.29 um`
increases maximum diffusive H2 flux ~10x, keeping diffusive capacity above biological
rates). Strictly refuses uncalibrated CPM lattice units (`:cpm_lattice_sites` without
a calibrated `lattice_pitch_um`).
"""
function evaluate_syntrophic_spatial_diffusion(;
        interspecies_distance_um::Real = 0.29,
        coordinate_unit::Symbol = :micrometers,
        lattice_pitch_um::Union{Nothing,Real} = nothing,
        D_H2_m2_s::Real = 4.5e-9,
        delta_C_H2_mol_m3::Real = 5.0e-4,
        aggregate_volume_fraction::Real = 0.15,
        specific_exchange_area_m_inv::Real = 1.5e5)
    if coordinate_unit == :cpm_lattice_sites
        (lattice_pitch_um !== nothing && isfinite(lattice_pitch_um) && lattice_pitch_um > 0.0) ||
            throw(ArgumentError(
                "SYN5 Spatial Calibration Violation: uncalibrated CPM lattice sites cannot be promoted to physical micrometer interspecies distance without calibrated lattice_pitch_um"
            ))
        d_um = Float64(interspecies_distance_um) * Float64(lattice_pitch_um)
    elseif coordinate_unit == :micrometers
        d_um = Float64(interspecies_distance_um)
    else
        throw(ArgumentError("coordinate_unit must be :micrometers or :cpm_lattice_sites, got :$(coordinate_unit)"))
    end

    (isfinite(d_um) && d_um > 0.0 && D_H2_m2_s > 0.0 && delta_C_H2_mol_m3 >= 0.0 &&
     0.0 < aggregate_volume_fraction <= 1.0 && specific_exchange_area_m_inv > 0.0) ||
        throw(ArgumentError("Invalid physical diffusion parameters in evaluate_syntrophic_spatial_diffusion"))

    d_m = d_um * 1.0e-6
    d_dispersed_m = 5.30e-6
    d_coagg_m = 0.29e-6

    # J_max [mol / (m^3 * s)] = D_H2 * a_spec * phi_agg * (delta_C / d)
    coeff = Float64(D_H2_m2_s) * Float64(specific_exchange_area_m_inv) * Float64(aggregate_volume_fraction) * Float64(delta_C_H2_mol_m3)
    j_max_mol_m3_s = coeff / d_m
    j_dispersed_mol_m3_s = coeff / d_dispersed_m
    j_coagg_mol_m3_s = coeff / d_coagg_m

    # Convert 1 mol/(m^3*s) = 60,000 nmol/(mL*min) to compare against Felchner-Zwirello 2013 (1.1 -> 10.3 nmol/(mL*min))
    j_max_nmol_mL_min = j_max_mol_m3_s * 60000.0
    flux_amplification = d_dispersed_m / d_m

    return (
        stage_id = :SYN5_spatial_aggregate_fickian_diffusion,
        coordinate_unit = coordinate_unit,
        calibrated_to_physical_meters = true,
        interspecies_distance_um = d_um,
        interspecies_distance_m = d_m,
        D_H2_m2_s = Float64(D_H2_m2_s),
        max_fickian_h2_flux_mol_m3_s = j_max_mol_m3_s,
        max_fickian_h2_flux_nmol_mL_min = j_max_nmol_mL_min,
        flux_at_5p30_um_mol_m3_s = j_dispersed_mol_m3_s,
        flux_at_0p29_um_mol_m3_s = j_coagg_mol_m3_s,
        fickian_flux_amplification_vs_5p30_um = flux_amplification,
        felchner_zwirello_2013_bounds_nmol_mL_min = (1.1, 10.3),
        spatial_coaggregation_increases_flux = (j_coagg_mol_m3_s > j_dispersed_mol_m3_s)
    )
end

"""
    export_syntrophic_eda_bundle_json(;
        output_path::AbstractString = joinpath(dirname(dirname(@__DIR__)), "matlab_fixtures", "syntrophy_demo_bundle.json"))

SYN6: Executes `SYN0..SYN5` and writes a schema-compliant `m7-field-bundle` v1.0.0 JSON
demonstration bundle (`syntrophy_demo_bundle.json`) for the read-only MATLAB R2026b
Syntrophy Reactor Explorer.
"""
function export_syntrophic_eda_bundle_json(;
        output_path::AbstractString = joinpath(dirname(dirname(@__DIR__)), "matlab_fixtures", "syntrophy_demo_bundle.json"))
    syn0 = evaluate_syntrophic_stoichiometric_ledger()
    syn0_butanol = evaluate_syntrophic_stoichiometric_ledger(
        substrate = :butyrate,
        donor_consumed_mol = 0.020,
        p_h2_ppm = 480.0,
        bicarbonate_limited = true,
        butanol_produced_mol = 0.004
    )
    syn1_37 = step_syntrophic_two_population_kinetics(temperature_C = 37.0, kinetic_model = :moser)
    syn1_65 = step_syntrophic_two_population_kinetics(temperature_C = 65.0, kinetic_model = :moser)
    syn1_70 = step_syntrophic_two_population_kinetics(temperature_C = 70.0, kinetic_model = :moser)
    syn2_phf = evaluate_syntrophic_gas_liquid_transfer(headspace_free_phf_hur = true)
    syn2_gas = evaluate_syntrophic_gas_liquid_transfer(headspace_free_phf_hur = false)
    syn3 = evaluate_syntrophic_thermodynamics(temperature_K = 310.15)
    syn4 = evaluate_syntrophic_identifiability_and_calibration()
    syn5_dispersed = evaluate_syntrophic_spatial_diffusion(interspecies_distance_um = 5.30)
    syn5_coagg = evaluate_syntrophic_spatial_diffusion(interspecies_distance_um = 0.29)

    temps_C = [27.0, 37.0, 50.0, 65.0, 70.0]
    monod_rates = Float64[]
    moser_rates = Float64[]
    qmax_vol_rates = Float64[]
    ks_vals = Float64[]
    for tc in temps_C
        st = step_syntrophic_two_population_kinetics(temperature_C = tc, kinetic_model = :moser, n_steps = 1)
        push!(monod_rates, st.monod_specific_rate_mol_kg_s * st.X_methanogen_final_kg_m3)
        push!(moser_rates, st.moser_specific_rate_mol_kg_s * st.X_methanogen_final_kg_m3)
        push!(qmax_vol_rates, st.q_max_H2_mol_kg_s * st.X_methanogen_final_kg_m3)
        push!(ks_vals, st.K_S_H2_mol_m3)
    end

    dists_um = [5.30, 2.50, 1.00, 0.50, 0.29]
    fick_fluxes = [evaluate_syntrophic_spatial_diffusion(interspecies_distance_um = d).max_fickian_h2_flux_mol_m3_s for d in dists_um]

    bundle = (
        schema = "m7-field-bundle",
        schema_version = "1.0.0",
        input_family = "syntrophic_bioenergetics_smbr",
        producer = (
            repo = "SmallModularBioReactors.jl",
            commit = "local_working_tree",
            solver = "SyntrophicBioenergetics_SYN0_to_SYN6"
        ),
        clock = (
            type = "physical_seconds",
            value = 1200.0,
            unit = "s"
        ),
        mesh = (
            type = "cartesian_3d",
            dims = [5, 1, 1],
            bounds_m = [[0.0, 0.05], [0.0, 0.01], [0.0, 0.01]],
            lattice_pitch_m = 0.01
        ),
        fields = (
            X_donor_kg_m3 = (
                unit = "kg/m^3",
                basis = "qpcr_partner_biomass",
                dtype = "float64",
                inline_values = fill(syn1_37.X_donor_final_kg_m3, 5)
            ),
            X_methanogen_kg_m3 = (
                unit = "kg/m^3",
                basis = "qpcr_partner_biomass",
                dtype = "float64",
                inline_values = fill(syn1_37.X_methanogen_final_kg_m3, 5)
            ),
            dissolved_H2_mol_m3 = (
                unit = "mol/m^3",
                basis = "pore_water",
                dtype = "float64",
                inline_values = ks_vals
            ),
            CH4_produced_mol_m3 = (
                unit = "mol/m^3",
                basis = "bulk_reactor",
                dtype = "float64",
                inline_values = fill(syn1_37.CH4_produced_mol_m3, 5)
            ),
            moser_H2_uptake_rate = (
                unit = "mol/(m^3*s)",
                basis = "phf_hur_headspace_free",
                dtype = "float64",
                inline_values = moser_rates
            ),
            monod_H2_uptake_rate = (
                unit = "mol/(m^3*s)",
                basis = "phf_hur_headspace_free",
                dtype = "float64",
                inline_values = monod_rates
            ),
            fickian_H2_diffusive_flux = (
                unit = "mol/(m^3*s)",
                basis = "coaggregate_calibrated_um",
                dtype = "float64",
                inline_values = fick_fluxes
            ),
            interspecies_distance_m = (
                unit = "m",
                basis = "fish_3d_calibrated",
                dtype = "float64",
                inline_values = dists_um .* 1.0e-6
            ),
            delta_r_G_donor_J_mol = (
                unit = "J/mol",
                basis = "transformed_thermodynamic",
                dtype = "float64",
                inline_values = fill(syn3.delta_r_G_donor_J_mol, 5)
            ),
            delta_r_G_methanogen_J_mol = (
                unit = "J/mol",
                basis = "transformed_thermodynamic",
                dtype = "float64",
                inline_values = fill(syn3.delta_r_G_methanogen_J_mol, 5)
            ),
            temperature_K = (
                unit = "K",
                basis = "reactor_thermal",
                dtype = "float64",
                inline_values = temps_C .+ 273.15
            ),
        ),
        evidence = (
            status = "Incomplete",
            calibration = "synthetic_demonstration_from_published_priors",
            evolutionary_track_verdict = "Incomplete",
            experimental_reactor_track_verdict = "Accepted_Model_Verification_Only"
        ),
        syntrophy_stages = (
            SYN0_stoichiometry = (
                propionate_e_rel_residual = syn0.electron_balance_rel_residual,
                propionate_e_closed = syn0.electron_balance_closed,
                butanol_shift_pH2_ppm = syn0_butanol.p_h2_ppm,
                butanol_shift_active = syn0_butanol.butanol_shift_active,
                butanol_produced_mol = syn0_butanol.butanol_produced_mol
            ),
            SYN1_phf_hur_kinetics = (
                temps_C = temps_C,
                monod_rates_mol_m3_s = monod_rates,
                moser_rates_mol_m3_s = moser_rates,
                qmax_volumetric_mol_m3_s = qmax_vol_rates,
                q_max_65C_mmol_gVSS_d = syn1_65.q_max_H2_mmol_gVSS_d,
                q_max_70C_mmol_gVSS_d = syn1_70.q_max_H2_mmol_gVSS_d,
                thermal_cutoff_70C_active = syn1_70.thermal_cutoff_at_70C,
                sse_monod = syn1_37.sse_monod,
                sse_moser = syn1_37.sse_moser,
                moser_preferred_over_monod = syn1_37.moser_preferred_over_monod
            ),
            SYN2_gas_liquid_ks_bias = (
                ks_true_mol_m3 = syn2_gas.K_S_true_mol_m3,
                ks_apparent_gas_supplied_mol_m3 = syn2_gas.K_S_apparent_gas_supplied_mol_m3,
                ks_bias_ratio_gas_supplied = syn2_gas.ks_bias_ratio,
                ks_bias_ratio_phf_hur = syn2_phf.ks_bias_ratio,
                oversaturation_ratio = syn2_gas.oversaturation_ratio
            ),
            SYN3_thermodynamics = (
                delta_r_G_donor_J_mol = syn3.delta_r_G_donor_J_mol,
                delta_r_G_methanogen_J_mol = syn3.delta_r_G_methanogen_J_mol,
                F_T_donor = syn3.thermodynamic_factor_F_T_donor,
                F_T_methanogen = syn3.thermodynamic_factor_F_T_methanogen,
                regulatory_silence_separate = true
            ),
            SYN4_identifiability_and_quarantine = (
                fisher_rank = Int64(syn4.fisher_rank),
                n_parameters = Int64(syn4.n_parameters),
                identifiable = syn4.identifiable,
                fisher_condition_number = syn4.fisher_condition_number,
                hydrogen_hypothesis_verdict = "Incomplete",
                asgard_heimdallarchaeia_verdict = "Incomplete",
                mcmc_ppc_rmse = syn4.mcmc_posterior.ppc_rmse
            ),
            SYN5_spatial_aggregate = (
                distances_um = dists_um,
                fickian_flux_mol_m3_s = fick_fluxes,
                flux_amplification_5p30_to_0p29 = syn5_coagg.fickian_flux_amplification_vs_5p30_um,
                flux_5p30_um = syn5_dispersed.max_fickian_h2_flux_mol_m3_s,
                flux_0p29_um = syn5_coagg.max_fickian_h2_flux_mol_m3_s
            )
        )
    )

    mkpath(dirname(output_path))
    json_txt = _to_json_str(bundle) * "\n"
    (!occursin('\u2014', json_txt) && !occursin('\u2013', json_txt)) ||
        throw(ArgumentError("Typography violation: em-dash or en-dash detected in syntrophy bundle JSON"))
    write(output_path, json_txt)
    return (
        bundle_path = String(output_path),
        bundle_sha256 = bytes2hex(sha256(Vector{UInt8}(json_txt))),
        bundle = bundle
    )
end
