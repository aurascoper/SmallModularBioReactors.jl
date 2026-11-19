# ==============================================================================
# src/transport/ClayBiofilmMetalInterface.jl
# M7.9 Clay-Biofilm-Metal Geochemical Interface (CBM-1):
# - Measurable 5-field primary vector X = [theta_mobile, aw, Fe2, Fe3, S_aq]
# - 5-pool sulfur ledger (SO4, S_aq, FeS/Cu2S, sorbed S, H2S(g))
# - 4-term clay density separation (SKB wet cutoff, SKB target, Grimsel dry, unresolved triad)
# - Discrete 0W..3W hydration & EC-TEM d001 shift + conventional vacuum TEM quarantine
# - Gates CBM-G1 (abiotic compaction/O2), CBM-G2 (screening vs bridging), CBM-G3 (DIET vs e-pili)
# - Fisher identifiability rank check & bead-calibrated FCM observation operator
# - Deterministic JSON field-bundle export for MATLAB R2026b 3D EDA
# ==============================================================================

const VALID_CLAY_HYDRATION_STATES = (
    :state_0W,
    :state_1W,
    :state_2W,
    :state_3W,
)

const VALID_CATION_POLYMER_MODES = (
    :charge_screening,
    :chelation,
    :polymer_bridging,
)

"""
    ClayBiofilmMetalState

M7.9 Clay-Biofilm-Metal Geochemical Interface state across `N` cells carrying the
five primary measurable fields `X = [theta_mobile, aw, Fe2, Fe3, S_aq]` alongside
EPS, exchangeable `Ca2+`/`Na+`, sulfide precipitates, magnetite, attached biomass,
and interfacial corrosion depth.
"""
struct ClayBiofilmMetalState
    cell_count::Int
    theta_mobile::Vector{Float64}
    aw::Vector{Float64}
    Fe2_mol_m3::Vector{Float64}
    Fe3_mol_m3::Vector{Float64}
    S_aq_mol_m3::Vector{Float64}
    EPS_kg_m3::Vector{Float64}
    Ca2_mol_m3::Vector{Float64}
    FeS_mineral_mol_m3::Vector{Float64}
    Fe3O4_magnetite_mol_m3::Vector{Float64}
    biomass_attached_kg_m3::Vector{Float64}
    corrosion_depth_m::Vector{Float64}

    function ClayBiofilmMetalState(;
            cell_count::Int = 5,
            theta_mobile::AbstractVector{<:Real} = [0.42, 0.38, 0.34, 0.30, 0.26],
            aw::AbstractVector{<:Real} = [0.992, 0.983, 0.972, 0.954, 0.935],
            Fe2_mol_m3::AbstractVector{<:Real} = [0.15, 0.28, 0.42, 0.55, 0.62],
            Fe3_mol_m3::AbstractVector{<:Real} = [2.40, 2.25, 2.10, 1.95, 1.85],
            S_aq_mol_m3::AbstractVector{<:Real} = [0.02, 0.05, 0.08, 0.04, 0.01],
            EPS_kg_m3::AbstractVector{<:Real} = [0.18, 0.22, 0.25, 0.19, 0.12],
            Ca2_mol_m3::AbstractVector{<:Real} = [1.50, 1.65, 1.80, 2.10, 2.40],
            FeS_mineral_mol_m3::AbstractVector{<:Real} = [0.05, 0.12, 0.18, 0.22, 0.15],
            Fe3O4_magnetite_mol_m3::AbstractVector{<:Real} = [0.02, 0.04, 0.06, 0.08, 0.10],
            biomass_attached_kg_m3::AbstractVector{<:Real} = [0.65, 0.72, 0.68, 0.45, 0.25],
            corrosion_depth_m::AbstractVector{<:Real} = [1.2e-6, 1.8e-6, 2.1e-6, 2.6e-6, 3.1e-6])
        cell_count >= 1 || throw(ArgumentError("cell_count must be >= 1"))
        vecs = (theta_mobile, aw, Fe2_mol_m3, Fe3_mol_m3, S_aq_mol_m3,
                EPS_kg_m3, Ca2_mol_m3, FeS_mineral_mol_m3, Fe3O4_magnetite_mol_m3,
                biomass_attached_kg_m3, corrosion_depth_m)
        all(v -> length(v) == cell_count, vecs) ||
            throw(ArgumentError("All ClayBiofilmMetalState field vectors must have length $(cell_count)"))
        all(x -> isfinite(x) && 0.0 < x < 1.0, theta_mobile) ||
            throw(ArgumentError("theta_mobile must lie strictly in (0, 1)"))
        all(x -> isfinite(x) && 0.0 < x <= 1.0, aw) ||
            throw(ArgumentError("water activity aw must lie in (0, 1]"))
        for v in (Fe2_mol_m3, Fe3_mol_m3, S_aq_mol_m3, EPS_kg_m3, Ca2_mol_m3,
                  FeS_mineral_mol_m3, Fe3O4_magnetite_mol_m3, biomass_attached_kg_m3, corrosion_depth_m)
            all(x -> isfinite(x) && x >= 0.0, v) ||
                throw(ArgumentError("Concentrations, mineral inventories, biomass, and corrosion depth must be finite and >= 0"))
        end
        new(
            cell_count,
            Float64.(theta_mobile),
            Float64.(aw),
            Float64.(Fe2_mol_m3),
            Float64.(Fe3_mol_m3),
            Float64.(S_aq_mol_m3),
            Float64.(EPS_kg_m3),
            Float64.(Ca2_mol_m3),
            Float64.(FeS_mineral_mol_m3),
            Float64.(Fe3O4_magnetite_mol_m3),
            Float64.(biomass_attached_kg_m3),
            Float64.(corrosion_depth_m)
        )
    end
end

"""
    evaluate_cbm_five_pool_sulfur_ledger(;
        initial_sulfate_aq_mol::Real = 1.0e-3,
        initial_sulfide_aq_mol::Real = 0.0,
        initial_precipitated_metal_sulfide_mol::Real = 0.0,
        initial_sorbed_sulfur_mol::Real = 0.0,
        initial_gas_h2s_mol::Real = 0.0,
        final_sulfate_aq_mol::Real = 6.2e-4,
        final_sulfide_aq_mol::Real = 5.0e-5,
        final_precipitated_metal_sulfide_mol::Real = 2.6e-4,
        final_sorbed_sulfur_mol::Real = 5.5e-5,
        final_gas_h2s_mol::Real = 1.5e-5,
        net_external_sulfur_input_mol::Real = 0.0,
        equate_fe_mineral_retention_with_system_removal::Bool = false,
        rel_tol::Real = 1.0e-6)

Audits the 5-pool sulfur mass balance distinguishing dissolved sulfate (`SO4^2-`),
dissolved sulfide (`S_aq`), precipitated metal sulfides (`FeS`/`Cu2S`), sorbed sulfur,
and gas-phase `H2S(g)`. Refuses equating iron-mineral sulfide retention with removal
from the system mass inventory.
"""
function evaluate_cbm_five_pool_sulfur_ledger(;
        initial_sulfate_aq_mol::Real = 1.0e-3,
        initial_sulfide_aq_mol::Real = 0.0,
        initial_precipitated_metal_sulfide_mol::Real = 0.0,
        initial_sorbed_sulfur_mol::Real = 0.0,
        initial_gas_h2s_mol::Real = 0.0,
        final_sulfate_aq_mol::Real = 6.2e-4,
        final_sulfide_aq_mol::Real = 5.0e-5,
        final_precipitated_metal_sulfide_mol::Real = 2.6e-4,
        final_sorbed_sulfur_mol::Real = 5.5e-5,
        final_gas_h2s_mol::Real = 1.5e-5,
        net_external_sulfur_input_mol::Real = 0.0,
        equate_fe_mineral_retention_with_system_removal::Bool = false,
        rel_tol::Real = 1.0e-6)
    equate_fe_mineral_retention_with_system_removal && throw(ArgumentError(
        "CBM Sulfur Accounting Violation: sulfide retention by iron minerals (FeS) is NOT equivalent to sulfur removal from the closed system mass balance"
    ))
    vals = (initial_sulfate_aq_mol, initial_sulfide_aq_mol, initial_precipitated_metal_sulfide_mol,
            initial_sorbed_sulfur_mol, initial_gas_h2s_mol, final_sulfate_aq_mol, final_sulfide_aq_mol,
            final_precipitated_metal_sulfide_mol, final_sorbed_sulfur_mol, final_gas_h2s_mol)
    all(v -> isfinite(v) && v >= 0.0, vals) ||
        throw(ArgumentError("All five sulfur pool molar inventories must be finite and non-negative"))
    (isfinite(net_external_sulfur_input_mol) && rel_tol > 0.0) ||
        throw(ArgumentError("Invalid net_external_sulfur_input_mol or rel_tol"))

    s_init = Float64(initial_sulfate_aq_mol) + Float64(initial_sulfide_aq_mol) +
             Float64(initial_precipitated_metal_sulfide_mol) + Float64(initial_sorbed_sulfur_mol) +
             Float64(initial_gas_h2s_mol)
    s_final = Float64(final_sulfate_aq_mol) + Float64(final_sulfide_aq_mol) +
              Float64(final_precipitated_metal_sulfide_mol) + Float64(final_sorbed_sulfur_mol) +
              Float64(final_gas_h2s_mol)
    s_expected = s_init + Float64(net_external_sulfur_input_mol)
    abs_res = abs(s_final - s_expected)
    rel_res = abs_res / max(1.0e-18, max(s_init, s_expected))

    sulfide_retained_fraction = (Float64(final_precipitated_metal_sulfide_mol) + Float64(final_sorbed_sulfur_mol)) /
                                max(1.0e-18, Float64(final_sulfide_aq_mol) + Float64(final_precipitated_metal_sulfide_mol) + Float64(final_sorbed_sulfur_mol) + Float64(final_gas_h2s_mol))

    return (
        pool_names = (:sulfate_aq, :sulfide_aq, :precipitated_metal_sulfide, :sorbed_sulfur, :gas_h2s),
        initial_total_S_mol = s_init,
        final_total_S_mol = s_final,
        expected_total_S_mol = s_expected,
        mass_balance_abs_residual_mol = abs_res,
        mass_balance_rel_residual = rel_res,
        sulfur_conserved = (rel_res <= Float64(rel_tol)),
        iron_mineral_retention_distinct_from_system_removal = true,
        retained_reduced_sulfur_fraction = sulfide_retained_fraction,
        final_pools_mol = (
            sulfate_aq = Float64(final_sulfate_aq_mol),
            sulfide_aq = Float64(final_sulfide_aq_mol),
            precipitated_metal_sulfide = Float64(final_precipitated_metal_sulfide_mol),
            sorbed_sulfur = Float64(final_sorbed_sulfur_mol),
            gas_h2s = Float64(final_gas_h2s_mol)
        )
    )
end

"""
    evaluate_clay_density_four_terms(; collapse_density_terms_to_single_threshold::Bool = false)

Enforces the four-term bentonite density distinction from `Clay, Microbes and Metal.pdf` p. 19:
1. `skb_tr1609_activity_cutoff_wet_kg_m3`: `1740..1880 kg/m^3` wet (`~1.2..1.4 g/cm^3` dry), measured.
2. `skb_design_target_wet_kg_m3`: `1950..2050 kg/m^3` wet design guidance (not a biological threshold).
3. `grimsel_field_measured_dry`: `a_w = 0.983` at `1.25 g/cm^3` dry and `0.954` at `1.50 g/cm^3` dry.
4. `unresolved_secondary_triad`: `(1.6 g/cm^3 dry, a_w 0.96, 2 MPa)` held in `:Incomplete`.
Throws `ArgumentError` if `collapse_density_terms_to_single_threshold == true`.
"""
function evaluate_clay_density_four_terms(; collapse_density_terms_to_single_threshold::Bool = false)
    collapse_density_terms_to_single_threshold && throw(ArgumentError(
        "CBM Density Terminology Violation (Clay PDF p. 19): wet density cut-off, SKB design target, Grimsel dry density, and unresolved secondary triad cannot be collapsed into one critical density parameter"
    ))

    return (
        four_terms_kept_separate = true,
        skb_tr1609_activity_cutoff_wet_kg_m3 = (
            range_wet_kg_m3 = (1740.0, 1880.0),
            approx_dry_g_cm3 = (1.2, 1.4),
            basis = :skb_tr1609_measured,
            is_biological_threshold = true,
            verdict = :Accepted
        ),
        skb_design_target_wet_kg_m3 = (
            range_wet_kg_m3 = (1950.0, 2050.0),
            basis = :engineering_design_guidance,
            is_biological_threshold = false,
            verdict = :Accepted
        ),
        grimsel_field_measured_dry = (
            point_1p25_dry_g_cm3 = (dry_g_cm3 = 1.25, aw = 0.983),
            point_1p50_dry_g_cm3 = (dry_g_cm3 = 1.50, aw = 0.954),
            basis = :grimsel_field_measured,
            verdict = :Accepted
        ),
        unresolved_secondary_triad = (
            dry_g_cm3 = 1.60,
            aw = 0.96,
            swelling_pressure_MPa = 2.0,
            basis = :secondary_citation_primary_text_unverified,
            verdict = :Incomplete
        )
    )
end

"""
    evaluate_clay_hydration_and_ectem(;
        hydration_state::Symbol = :state_2W,
        cation::Symbol = :Ca2_plus,
        relative_humidity_pct::Real = 40.0,
        fe3_bioreduced_fraction::Real = 0.0,
        imaging_mode::Symbol = :ec_tem_hydrated,
        promote_vacuum_tem_d001_as_hydrated::Bool = false)

Models discrete montmorillonite hydration (`0W..3W`, Ferrage et al. 2005a/b & Greenspan 1977)
and Environmental Cell TEM (`EC-TEM`, Dong et al. 2003) `d001` contraction (`1.5 -> 1.1 nm`)
upon structural Fe(III) reduction. Conventional high-vacuum TEM dehydrates smectite and is
quarantined at `:Incomplete`.
"""
function evaluate_clay_hydration_and_ectem(;
        hydration_state::Symbol = :state_2W,
        cation::Symbol = :Ca2_plus,
        relative_humidity_pct::Real = 40.0,
        fe3_bioreduced_fraction::Real = 0.0,
        imaging_mode::Symbol = :ec_tem_hydrated,
        promote_vacuum_tem_d001_as_hydrated::Bool = false)
    promote_vacuum_tem_d001_as_hydrated && throw(ArgumentError(
        "CBM Vacuum TEM Quarantine Violation (Dong et al. 2003): conventional high-vacuum TEM dehydrates smectite and collapses d001; cannot be promoted as hydrated d001 evidence"
    ))
    hydration_state in VALID_CLAY_HYDRATION_STATES ||
        throw(ArgumentError("hydration_state must be in $(VALID_CLAY_HYDRATION_STATES)"))
    cation in (:Na_plus, :Ca2_plus) ||
        throw(ArgumentError("cation must be :Na_plus or :Ca2_plus"))
    imaging_mode in (:ec_tem_hydrated, :xrd_rh_controlled, :conventional_vacuum_tem) ||
        throw(ArgumentError("imaging_mode must be :ec_tem_hydrated, :xrd_rh_controlled, or :conventional_vacuum_tem"))
    (0.0 <= relative_humidity_pct <= 100.0 && 0.0 <= fe3_bioreduced_fraction <= 1.0) ||
        throw(ArgumentError("relative_humidity_pct must lie in [0, 100] and fe3_bioreduced_fraction in [0, 1]"))

    baseline_d001_A = if hydration_state == :state_0W
        10.0
    elseif hydration_state == :state_1W
        12.6
    elseif hydration_state == :state_2W
        15.25
    else
        (cation == :Na_plus) ? 18.75 : 15.3 # osmotic 3W is Na-specific
    end

    # Dong et al. 2003 EC-TEM: unreduced 1.5 +/- 0.1 nm -> bioreduced 1.1 +/- 0.1 nm (at ~26% Fe(III) reduction)
    red_scale = clamp(Float64(fe3_bioreduced_fraction) / 0.26, 0.0, 1.0)
    baseline_d001_nm = baseline_d001_A / 10.0
    ectem_d001_nm = hydration_state == :state_2W ?
                    (1.5 - 0.4 * red_scale) :
                    (baseline_d001_nm * (1.0 - (0.4 / 1.5) * red_scale))

    if imaging_mode == :conventional_vacuum_tem
        d001_observed_nm = 1.00 # dehydrated vacuum collapse
        verdict = :Incomplete
        vacuum_quarantined = true
    else
        d001_observed_nm = (fe3_bioreduced_fraction > 0.0 && imaging_mode == :ec_tem_hydrated) ?
                           ectem_d001_nm : baseline_d001_nm
        verdict = :Accepted
        vacuum_quarantined = false
    end

    osmotic_swelling_possible = (cation == :Na_plus) && (hydration_state == :state_3W)

    return (
        hydration_state = hydration_state,
        cation = cation,
        relative_humidity_pct = Float64(relative_humidity_pct),
        baseline_d001_A = baseline_d001_A,
        fe3_bioreduced_fraction = Float64(fe3_bioreduced_fraction),
        ectem_d001_nm = ectem_d001_nm,
        observed_d001_nm = d001_observed_nm,
        observed_d001_m = d001_observed_nm * 1.0e-9,
        imaging_mode = imaging_mode,
        vacuum_tem_quarantined = vacuum_quarantined,
        osmotic_swelling_possible = osmotic_swelling_possible,
        verdict = verdict
    )
end

"""
    evaluate_cbm_g1_abiotic_compaction_gate(;
        dry_density_kg_m3::Real = 1600.0,
        abiotic_sterile_rate_um_yr::Real = 2.4,
        biotic_rate_um_yr::Real = 1.8,
        detection_limit_um_yr::Real = 0.25,
        has_sterile_control::Bool = true,
        has_killed_cell_control::Bool = true,
        o2_excluded_anoxic::Bool = true,
        assume_higher_biomass_or_lower_density_means_higher_corrosion::Bool = false)

CBM-G1: Evaluates abiotic vs biotic corrosion under bentonite compaction (Cerna et al. 2026:
sterile/low-biomass systems at high compaction can exhibit higher corrosion rates than
biotic treatments; O2 ingress must be excluded before attributing corrosion to MIC).
"""
function evaluate_cbm_g1_abiotic_compaction_gate(;
        dry_density_kg_m3::Real = 1600.0,
        abiotic_sterile_rate_um_yr::Real = 2.4,
        biotic_rate_um_yr::Real = 1.8,
        detection_limit_um_yr::Real = 0.25,
        has_sterile_control::Bool = true,
        has_killed_cell_control::Bool = true,
        o2_excluded_anoxic::Bool = true,
        assume_higher_biomass_or_lower_density_means_higher_corrosion::Bool = false)
    assume_higher_biomass_or_lower_density_means_higher_corrosion && throw(ArgumentError(
        "CBM-G1 Violation (Cerna et al. 2026): higher microbial abundance or lower compaction cannot be assumed a priori to yield higher corrosion than sterile/compacted controls"
    ))
    (dry_density_kg_m3 > 0.0 && abiotic_sterile_rate_um_yr >= 0.0 && biotic_rate_um_yr >= 0.0 && detection_limit_um_yr >= 0.0) ||
        throw(ArgumentError("Density, corrosion rates, and detection limit must be valid non-negative numbers"))

    controls_ok = has_sterile_control && has_killed_cell_control && o2_excluded_anoxic
    abiotic = Float64(abiotic_sterile_rate_um_yr)
    biotic = Float64(biotic_rate_um_yr)
    lod = Float64(detection_limit_um_yr)
    excess_biotic = biotic - abiotic
    below_detection_limit = (abiotic < lod) && (biotic < lod)
    abiotic_dominated = !below_detection_limit && (abiotic >= biotic)

    verdict = if !controls_ok || below_detection_limit
        :Incomplete # Uninformative: missing controls or both rates below gravimetric detection limit
    elseif abiotic_dominated
        :Rejected # Biotic-MIC dominance hypothesis rejected when sterile >= biotic above LOD
    elseif excess_biotic > lod
        :Accepted
    else
        :Incomplete # Excess below detection limit / uninformative
    end

    return (
        gate_id = :CBM_G1_abiotic_compaction_and_O2_gate,
        dry_density_kg_m3 = Float64(dry_density_kg_m3),
        abiotic_sterile_rate_um_yr = abiotic,
        biotic_rate_um_yr = biotic,
        excess_biotic_rate_um_yr = excess_biotic,
        below_detection_limit = below_detection_limit,
        abiotic_dominated = abiotic_dominated,
        controls_satisfied = controls_ok,
        biotic_mic_attribution_verdict = verdict
    )
end

"""
    evaluate_cbm_g2_screening_vs_bridging_gate(;
        proposed_mode::Symbol = :polymer_bridging,
        matched_monovalent_control::Bool = true,
        require_empirical_rheology::Bool = true,
        measured_rheology_or_retention_shift::Bool = false,
        zeta_potential_shifted::Bool = true,
        infer_bridging_from_zeta_potential_only::Bool = false,
        equate_prescribed_crosslink_fn_with_validation::Bool = false)

CBM-G2: Separates charge screening, IUPAC chelation (`C01012`), and true polymer bridging.
Bridging requires BOTH a matched ionic-strength/pH monovalent control AND empirical
rheology/retention evidence (`require_empirical_rheology = true`, `measured_rheology_or_retention_shift = true`).
"""
function evaluate_cbm_g2_screening_vs_bridging_gate(;
        proposed_mode::Symbol = :polymer_bridging,
        matched_monovalent_control::Bool = true,
        require_empirical_rheology::Bool = true,
        measured_rheology_or_retention_shift::Bool = false,
        zeta_potential_shifted::Bool = true,
        infer_bridging_from_zeta_potential_only::Bool = false,
        equate_prescribed_crosslink_fn_with_validation::Bool = false)
    infer_bridging_from_zeta_potential_only && throw(ArgumentError(
        "CBM-G2 Violation (Clay PDF p. 8): surface/zeta potential shift with ionic strength demonstrates charge screening, not EPS cross-linking"
    ))
    equate_prescribed_crosslink_fn_with_validation && throw(ArgumentError(
        "CBM-G2 Violation: prescribed PolymerTransport.jl Ca-crosslink function is a model assumption, not empirical clay-EPS bridging validation"
    ))
    proposed_mode in VALID_CATION_POLYMER_MODES ||
        throw(ArgumentError("proposed_mode must be in $(VALID_CATION_POLYMER_MODES)"))

    bridging_verified = (proposed_mode == :polymer_bridging) &&
                        matched_monovalent_control &&
                        require_empirical_rheology &&
                        measured_rheology_or_retention_shift

    verdict = if proposed_mode == :polymer_bridging
        bridging_verified ? :Accepted : :Incomplete
    elseif proposed_mode == :charge_screening
        zeta_potential_shifted ? :Accepted : :Incomplete
    else
        :Incomplete
    end

    return (
        gate_id = :CBM_G2_screening_chelation_bridging_gate,
        proposed_mode = proposed_mode,
        matched_monovalent_control = matched_monovalent_control,
        require_empirical_rheology = require_empirical_rheology,
        measured_rheology_or_retention_shift = measured_rheology_or_retention_shift,
        bridging_verdict = (proposed_mode == :polymer_bridging ? verdict : :Incomplete),
        verdict = verdict
    )
end

"""
    evaluate_cbm_g3_diet_vs_epili_gate(;
        magnetite_present::Bool = true,
        e_pili_present_and_verified::Bool = false,
        soluble_shuttles_excluded::Bool = false,
        direct_flux_perturbation_verified::Bool = false,
        ch4_stimulated_in_epili_deficient_control::Bool = false,
        claim_diet_from_fe3o4_presence_alone::Bool = false)

CBM-G3: Enforces Wang et al. (2018) & Cruz Viggi et al. (2014): `Fe3O4` detection alone
does not establish DIET (`:Incomplete`), and magnetite stimulation cannot substitute for
electrically conductive pili (e-pili) unless verified against an e-pili-deficient control
and soluble shuttle exclusion.
"""
function evaluate_cbm_g3_diet_vs_epili_gate(;
        magnetite_present::Bool = true,
        e_pili_present_and_verified::Bool = false,
        soluble_shuttles_excluded::Bool = false,
        direct_flux_perturbation_verified::Bool = false,
        ch4_stimulated_in_epili_deficient_control::Bool = false,
        claim_diet_from_fe3o4_presence_alone::Bool = false)
    claim_diet_from_fe3o4_presence_alone && throw(ArgumentError(
        "CBM-G3 DIET Violation (Wang et al. 2018; Clay PDF): Fe3O4 magnetite detection alone does not establish Direct Interspecies Electron Transfer"
    ))

    verdict = if magnetite_present && e_pili_present_and_verified && soluble_shuttles_excluded && direct_flux_perturbation_verified
        :Accepted
    elseif magnetite_present && !e_pili_present_and_verified && !ch4_stimulated_in_epili_deficient_control && direct_flux_perturbation_verified
        :Rejected # Wang 2018 pattern: magnetite fails to stimulate CH4 without e-pili
    else
        :Incomplete
    end

    return (
        gate_id = :CBM_G3_magnetite_diet_vs_epili_gate,
        magnetite_present = magnetite_present,
        e_pili_present_and_verified = e_pili_present_and_verified,
        soluble_shuttles_excluded = soluble_shuttles_excluded,
        direct_flux_perturbation_verified = direct_flux_perturbation_verified,
        diet_causal_verdict = verdict
    )
end

"""
    evaluate_cbm_identifiability_and_fcm_operator(;
        measure_aw_and_d001::Bool = true,
        measure_fe_speciation::Bool = true,
        measure_sulfide_and_avs::Bool = true,
        measure_bead_calibrated_fcm::Bool = true,
        cell_events::Integer = 4500,
        bead_events::Integer = 500,
        bead_concentration_per_mL::Real = 1.0e5,
        dry_mass_per_cell_kg::Real = 2.8e-16,
        has_clay_debris_and_dead_controls::Bool = true,
        use_uncalibrated_raw_fcm_events::Bool = false,
        claim_spatial_position_from_fcm::Bool = false)

Evaluates the M7.9 parameter-identifiability study (hydration limitation vs donor
limitation vs Fe(III) reduction vs sulfide retention) via Fisher Information Matrix rank,
and applies the bead-calibrated Flow Cytometry (`FCM`) observation operator while enforcing
spatial-position, EPS cross-link, and equifinality boundaries.
"""
function evaluate_cbm_identifiability_and_fcm_operator(;
        measure_aw_and_d001::Bool = true,
        measure_fe_speciation::Bool = true,
        measure_sulfide_and_avs::Bool = true,
        measure_bead_calibrated_fcm::Bool = true,
        cell_events::Integer = 4500,
        bead_events::Integer = 500,
        bead_concentration_per_mL::Real = 1.0e5,
        dry_mass_per_cell_kg::Real = 2.8e-16,
        has_clay_debris_and_dead_controls::Bool = true,
        use_uncalibrated_raw_fcm_events::Bool = false,
        claim_spatial_position_from_fcm::Bool = false)
    use_uncalibrated_raw_fcm_events && throw(ArgumentError(
        "CBM FCM Violation: raw uncalibrated FCM event counts cannot serve as absolute suspended cell concentration without counting-bead calibration"
    ))
    claim_spatial_position_from_fcm && throw(ArgumentError(
        "CBM FCM Scope Violation: flow cytometry measures suspended/detached phenotypes and cannot resolve spatial position inside clay pores or EPS cross-link density"
    ))
    (cell_events >= 0 && bead_events >= 1 && bead_concentration_per_mL > 0.0 && dry_mass_per_cell_kg > 0.0) ||
        throw(ArgumentError("Invalid FCM event/bead parameters"))

    # Parameter vector theta = [k_hydration_aw, k_donor_lim, k_fe3_red, k_sulfide_ret] (4 params)
    rows = Vector{Vector{Float64}}()
    # Bulk effluent/sulfate-loss row alone confounds hydration & donor limitation without Fe(III)/AVS resolution
    push!(rows, [0.50, 0.50, 0.00, 0.00])
    if measure_aw_and_d001
        push!(rows, [1.00, 0.05, 0.00, 0.00])
    end
    if measure_fe_speciation
        push!(rows, [0.00, 0.15, 1.00, 0.10])
    end
    if measure_sulfide_and_avs
        push!(rows, [0.00, 0.10, 0.15, 1.00])
    end
    if measure_bead_calibrated_fcm
        if measure_aw_and_d001
            push!(rows, [0.15, 0.90, 0.00, 0.00])
        else
            push!(rows, [0.50, 0.50, 0.00, 0.00])
        end
    end

    J = reduce(vcat, [r' for r in rows])
    F = J' * J
    sv = svdvals(F)
    tol_rank = 1.0e-8 * max(1.0, maximum(sv))
    rank_F = count(s -> s > tol_rank, sv)
    n_params = 4
    identifiable = (rank_F == n_params)
    cond_num = minimum(sv) > 1.0e-15 ? (maximum(sv) / minimum(sv)) : Inf

    # Bead-calibrated FCM observation operator:
    # C_cells [cells/mL] = (cell_events / bead_events) * bead_concentration_per_mL
    # Conversion to SI kg/m^3: C_cells * 1e6 [cells/m^3] * dry_mass_per_cell_kg
    cells_per_mL = (Float64(cell_events) / Float64(bead_events)) * Float64(bead_concentration_per_mL)
    cells_per_m3 = cells_per_mL * 1.0e6
    fcm_biomass_kg_m3 = cells_per_m3 * Float64(dry_mass_per_cell_kg)
    fcm_verdict = (measure_bead_calibrated_fcm && has_clay_debris_and_dead_controls) ? :Accepted : :Incomplete

    return (
        n_parameters = n_params,
        fisher_rank = rank_F,
        identifiable = identifiable,
        fisher_eigenvalues = sv,
        fisher_condition_number = cond_num,
        fcm_observation_operator = (
            cells_per_mL = cells_per_mL,
            cells_per_m3 = cells_per_m3,
            suspended_biomass_kg_m3 = fcm_biomass_kg_m3,
            bead_calibrated = measure_bead_calibrated_fcm,
            clay_debris_and_dead_controls_present = has_clay_debris_and_dead_controls,
            resolves_clay_pore_spatial_position = false,
            measures_eps_crosslink_density = false,
            equifinality_caveat_applies = true,
            verdict = fcm_verdict
        )
    )
end

"""
    export_cbm_eda_bundle_json(;
        output_path::AbstractString = joinpath(dirname(dirname(@__DIR__)), "matlab_fixtures", "cbm_m79_demo_bundle.json"))

Generates and writes `cbm_m79_demo_bundle.json` (`m7-field-bundle` v1.0.0) for the read-only
MATLAB R2026b Clay-Biofilm-Metal Explorer (`M7.9-CBM`).
"""
function export_cbm_eda_bundle_json(;
        output_path::AbstractString = joinpath(dirname(dirname(@__DIR__)), "matlab_fixtures", "cbm_m79_demo_bundle.json"))
    st = ClayBiofilmMetalState()
    s_ledger = evaluate_cbm_five_pool_sulfur_ledger()
    dens4 = evaluate_clay_density_four_terms()
    hyd_2w = evaluate_clay_hydration_and_ectem(hydration_state = :state_2W, cation = :Ca2_plus, imaging_mode = :xrd_rh_controlled)
    hyd_ectem = evaluate_clay_hydration_and_ectem(hydration_state = :state_2W, fe3_bioreduced_fraction = 0.26, imaging_mode = :ec_tem_hydrated)
    hyd_vac = evaluate_clay_hydration_and_ectem(hydration_state = :state_2W, imaging_mode = :conventional_vacuum_tem)
    g1 = evaluate_cbm_g1_abiotic_compaction_gate()
    g2_inc = evaluate_cbm_g2_screening_vs_bridging_gate(proposed_mode = :polymer_bridging, require_empirical_rheology = false, measured_rheology_or_retention_shift = false)
    g2_acc = evaluate_cbm_g2_screening_vs_bridging_gate(proposed_mode = :polymer_bridging, require_empirical_rheology = true, matched_monovalent_control = true, measured_rheology_or_retention_shift = true)
    g3_inc = evaluate_cbm_g3_diet_vs_epili_gate(magnetite_present = true, e_pili_present_and_verified = false)
    g3_rej = evaluate_cbm_g3_diet_vs_epili_gate(magnetite_present = true, e_pili_present_and_verified = false, direct_flux_perturbation_verified = true)
    ident_full = evaluate_cbm_identifiability_and_fcm_operator()
    ident_drop = evaluate_cbm_identifiability_and_fcm_operator(measure_aw_and_d001 = false, measure_sulfide_and_avs = false)

    d001_states_m = [
        evaluate_clay_hydration_and_ectem(hydration_state = :state_0W).observed_d001_m,
        evaluate_clay_hydration_and_ectem(hydration_state = :state_1W, cation = :Na_plus).observed_d001_m,
        hyd_2w.observed_d001_m,
        evaluate_clay_hydration_and_ectem(hydration_state = :state_3W, cation = :Na_plus).observed_d001_m,
        hyd_ectem.observed_d001_m,
    ]

    bundle = (
        schema = "m7-field-bundle",
        schema_version = "1.0.0",
        input_family = "clay_biofilm_metal_m79",
        producer = (
            repo = "SmallModularBioReactors.jl",
            commit = "local_working_tree",
            solver = "ClayBiofilmMetalInterface_M79_CBM"
        ),
        clock = (
            type = "physical_seconds",
            value = 1814400.0,
            unit = "s"
        ),
        mesh = (
            type = "cartesian_3d",
            dims = [5, 1, 1],
            bounds_m = [[0.0, 0.05], [0.0, 0.01], [0.0, 0.01]],
            lattice_pitch_m = 0.01
        ),
        fields = (
            theta_mobile = (
                unit = "1",
                basis = "mobile_pore_water",
                dtype = "float64",
                inline_values = st.theta_mobile
            ),
            aw = (
                unit = "1",
                basis = "clay_water_activity",
                dtype = "float64",
                inline_values = st.aw
            ),
            Fe2_mol_m3 = (
                unit = "mol/m^3",
                basis = "pore_and_exchangeable",
                dtype = "float64",
                inline_values = st.Fe2_mol_m3
            ),
            Fe3_mol_m3 = (
                unit = "mol/m^3",
                basis = "structural_and_reactive_iron",
                dtype = "float64",
                inline_values = st.Fe3_mol_m3
            ),
            S_aq_mol_m3 = (
                unit = "mol/m^3",
                basis = "dissolved_sulfide",
                dtype = "float64",
                inline_values = st.S_aq_mol_m3
            ),
            EPS_kg_m3 = (
                unit = "kg/m^3",
                basis = "biofilm_matrix",
                dtype = "float64",
                inline_values = st.EPS_kg_m3
            ),
            Ca2_mol_m3 = (
                unit = "mol/m^3",
                basis = "interfacial_cation",
                dtype = "float64",
                inline_values = st.Ca2_mol_m3
            ),
            FeS_mineral_mol_m3 = (
                unit = "mol/m^3",
                basis = "precipitated_sulfide_phase",
                dtype = "float64",
                inline_values = st.FeS_mineral_mol_m3
            ),
            d001_m = (
                unit = "m",
                basis = "clay_basal_spacing_0W_to_3W_and_ectem",
                dtype = "float64",
                inline_values = d001_states_m
            ),
            corrosion_depth_m = (
                unit = "m",
                basis = "interfacial_metal_loss",
                dtype = "float64",
                inline_values = st.corrosion_depth_m
            ),
        ),
        evidence = (
            status = "Incomplete",
            calibration = "synthetic_demonstration_awaiting_bsl1_pilot",
            note = "Read-only M7.9 demonstration; not an experimentally validated reactive-transport prediction"
        ),
        cbm_m79_stages = (
            five_pool_sulfur_ledger = (
                pools_mol = [
                    s_ledger.final_pools_mol.sulfate_aq,
                    s_ledger.final_pools_mol.sulfide_aq,
                    s_ledger.final_pools_mol.precipitated_metal_sulfide,
                    s_ledger.final_pools_mol.sorbed_sulfur,
                    s_ledger.final_pools_mol.gas_h2s
                ],
                pool_labels = ["SO4_aq", "S_aq", "FeS_Cu2S", "Sorbed_S", "H2S_gas"],
                sulfur_conserved = s_ledger.sulfur_conserved,
                rel_residual = s_ledger.mass_balance_rel_residual
            ),
            four_term_clay_density = (
                skb_cutoff_wet_kg_m3 = [1740.0, 1880.0],
                skb_design_target_wet_kg_m3 = [1950.0, 2050.0],
                grimsel_dry_g_cm3 = [1.25, 1.50],
                grimsel_aw = [0.983, 0.954],
                unresolved_triad_verdict = String(dens4.unresolved_secondary_triad.verdict)
            ),
            hydration_0W_to_3W_and_ectem = (
                labels = ["0W_dry", "1W_Na", "2W_Ca", "3W_Na", "EC_TEM_reduced"],
                d001_nm = d001_states_m .* 1.0e9,
                vacuum_tem_verdict = String(hyd_vac.verdict),
                ec_tem_verdict = String(hyd_ectem.verdict)
            ),
            discriminating_gates = (
                CBM_G1_abiotic_dominated = g1.abiotic_dominated,
                CBM_G1_biotic_mic_verdict = String(g1.biotic_mic_attribution_verdict),
                CBM_G1_sterile_vs_biotic_um_yr = [g1.abiotic_sterile_rate_um_yr, g1.biotic_rate_um_yr],
                CBM_G2_without_rheology_verdict = String(g2_inc.verdict),
                CBM_G2_with_rheology_and_monovalent_verdict = String(g2_acc.verdict),
                CBM_G3_magnetite_only_verdict = String(g3_inc.diet_causal_verdict),
                CBM_G3_epili_absent_refutation_verdict = String(g3_rej.diet_causal_verdict)
            ),
            identifiability_and_fcm = (
                full_fisher_rank = Int64(ident_full.fisher_rank),
                underdetermined_fisher_rank = Int64(ident_drop.fisher_rank),
                n_parameters = Int64(ident_full.n_parameters),
                full_singular_values = ident_full.fisher_eigenvalues,
                fcm_cells_per_mL = ident_full.fcm_observation_operator.cells_per_mL,
                fcm_biomass_kg_m3 = ident_full.fcm_observation_operator.suspended_biomass_kg_m3,
                fcm_resolves_clay_pore_position = ident_full.fcm_observation_operator.resolves_clay_pore_spatial_position,
                fcm_verdict = String(ident_full.fcm_observation_operator.verdict)
            )
        )
    )

    mkpath(dirname(output_path))
    json_txt = _to_json_str(bundle) * "\n"
    (!occursin('\u2014', json_txt) && !occursin('\u2013', json_txt)) ||
        throw(ArgumentError("Typography violation: em-dash or en-dash detected in CBM M7.9 bundle JSON"))
    write(output_path, json_txt)
    return (
        bundle_path = String(output_path),
        bundle_sha256 = bytes2hex(sha256(Vector{UInt8}(json_txt))),
        bundle = bundle
    )
end
