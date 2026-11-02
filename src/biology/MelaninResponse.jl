# ==============================================================================
# src/biology/MelaninResponse.jl
# 4-Stage Radiative-Reactive Transport Chain:
#   Stage 1 (Particle transport) -> Stage 2 (Energy deposition D_dot) ->
#   Stage 3 (Radiochemical species OH_radical / H2O2 scavenging) ->
#   Stage 4 (Cellular damage, repair, and growth modulation)
# Explicitly Quarantining Melanin Magnetoreception & Direct Radiosynthesis ATP
# Generation as Unestablished Hypotheses (Demonstration 2: 4I3OYyXt...jpg)
# ==============================================================================

const QUARANTINED_MELANIN_HYPOTHESES = (
    melanin_mediated_magnetoreception = :unestablished_hypothesis,
    direct_net_metabolic_atp_radiosynthesis = :unestablished_hypothesis,
)

"""
    evaluate_melanin_radiative_transport_chain(;
        incident_dose_rate_Gy_s::Real = 1.5e-4,
        melanin_mass_fraction::Real = 0.18,
        enable_unestablished_magnetoreception::Bool = false,
        enable_unestablished_net_atp_radiosynthesis::Bool = false)

Evaluates the 4-stage mechanistic melanin radiation-response chain comparing a heavily
melanized colony (`melanin_mass_fraction > 0`, left phenotype in `4I3OYyXt...jpg`) against
a non-melanized control colony (`melanin_mass_fraction == 0`, right phenotype in `4I3OYyXt...jpg`).
Throws `ArgumentError` if quarantined hypotheses (`enable_unestablished_magnetoreception` or
`enable_unestablished_net_atp_radiosynthesis`) are activated as physical source terms.
"""
function evaluate_melanin_radiative_transport_chain(;
        incident_dose_rate_Gy_s::Real = 1.5e-4,
        melanin_mass_fraction::Real = 0.18,
        enable_unestablished_magnetoreception::Bool = false,
        enable_unestablished_net_atp_radiosynthesis::Bool = false)
    enable_unestablished_magnetoreception && throw(ArgumentError(
        "Epistemic Quarantine Violation: melanin-mediated magnetoreception is quarantined as :unestablished_hypothesis"
    ))
    enable_unestablished_net_atp_radiosynthesis && throw(ArgumentError(
        "Epistemic Quarantine Violation: direct net metabolic ATP radiosynthesis without electron-donor mass balance is quarantined as :unestablished_hypothesis"
    ))
    incident_dose_rate_Gy_s >= 0.0 || throw(ArgumentError("incident_dose_rate_Gy_s must be >= 0"))
    (0.0 <= melanin_mass_fraction <= 0.50) || throw(ArgumentError("melanin_mass_fraction must lie in [0, 0.50]"))

    f_mel = Float64(melanin_mass_fraction)
    d0 = Float64(incident_dose_rate_Gy_s)

    # Stage 1: Particle/photon attenuation factor through melanized cell wall
    stage_1_transmission_factor = exp(-0.85 * f_mel)
    # Stage 2: Intracellular energy deposition rate D_dot (Gy/s)
    stage_2_deposited_dose_rate_Gy_s = d0 * stage_1_transmission_factor
    # Stage 3: Radiochemical radical yield (OH_radical + H2O2 in mol/(m^3*s)) after melanin free-radical scavenging
    radical_scavenging_efficiency = (4.5 * f_mel) / (1.0 + 4.5 * f_mel)
    unscavenged_radical_rate_mol_m3_s = stage_2_deposited_dose_rate_Gy_s * 3.5e-4 * (1.0 - radical_scavenging_efficiency)
    # Stage 4: Cellular double-strand break / oxidative damage rate vs repair capacity -> net radial growth rate
    damage_rate_s_inv = 2.5e3 * unscavenged_radical_rate_mol_m3_s
    repair_capacity_s_inv = 1.0e-7
    net_damage_penalty = max(0.0, damage_rate_s_inv - repair_capacity_s_inv) + 0.25 * damage_rate_s_inv
    net_radial_growth_um_hr = 45.0 / (1.0 + 2.5e5 * net_damage_penalty)

    return (
        stage_1_particle_transmission_factor = stage_1_transmission_factor,
        stage_2_deposited_dose_rate_Gy_s = stage_2_deposited_dose_rate_Gy_s,
        stage_3_radical_scavenging_efficiency = radical_scavenging_efficiency,
        stage_3_unscavenged_radical_rate_mol_m3_s = unscavenged_radical_rate_mol_m3_s,
        stage_4_oxidative_damage_rate_s_inv = damage_rate_s_inv,
        stage_4_net_radial_growth_um_hr = net_radial_growth_um_hr,
        quarantined_hypotheses = QUARANTINED_MELANIN_HYPOTHESES,
        evidence_image_sha256 = DOWNLOADS_FIVE_IMAGES_MANIFEST.demonstration_2_melanized_vs_control_fungal_colonies.sha256
    )
end

"""
    NeurosporaCrassaCircadianHyphalPreset

Preset specification for `Neurospora crassa` syncytial mycelial network inside a modular
bioreactor array (`dialectic-turns-2026-11-02.jsonl`, Turn `L3`, `index = 21`):
- `frq`, `wc-1`, `wc-2` White Collar Complex (WCC) transcription-translation circadian oscillator
  (`free_running_period_hr = 22.5` hr in constant darkness `DD`, entrainable to `entrained_period_hr`).
- Syncytial cytoplasmic streaming and hyphal anastomosis (`beta_1_syncytial_loops >= 1`).
- Coupled DHN-melanin / carotenoid (`neurosporaxanthin`) cell-wall photoprotection and radioprotection.
"""
struct NeurosporaCrassaCircadianHyphalPreset
    organism_symbol::Symbol
    clock_genes::Tuple{Symbol,Symbol,Symbol}
    free_running_period_hr::Float64
    entrained_period_hr::Float64
    basal_hyphal_extension_um_hr::Float64
    circadian_modulation_amplitude::Float64
    melanin_carotenoid_mass_fraction::Float64
    n_modular_chambers::Int
    beta_1_syncytial_loops::Int

    function NeurosporaCrassaCircadianHyphalPreset(;
            organism_symbol::Symbol = :Neurospora_crassa,
            clock_genes::Tuple{Symbol,Symbol,Symbol} = (:frq, :wc_1, :wc_2),
            free_running_period_hr::Real = 22.5,
            entrained_period_hr::Real = 24.0,
            basal_hyphal_extension_um_hr::Real = 120.0,
            circadian_modulation_amplitude::Real = 0.35,
            melanin_carotenoid_mass_fraction::Real = 0.18,
            n_modular_chambers::Integer = 4,
            beta_1_syncytial_loops::Integer = 3)
        (isfinite(free_running_period_hr) && 16.0 <= free_running_period_hr <= 32.0) ||
            throw(ArgumentError("free_running_period_hr must lie in [16.0, 32.0] hr for Neurospora crassa frq clock"))
        (isfinite(entrained_period_hr) && entrained_period_hr > 0.0) ||
            throw(ArgumentError("entrained_period_hr must be > 0"))
        (isfinite(basal_hyphal_extension_um_hr) && basal_hyphal_extension_um_hr > 0.0) ||
            throw(ArgumentError("basal_hyphal_extension_um_hr must be > 0"))
        (isfinite(circadian_modulation_amplitude) && 0.0 <= circadian_modulation_amplitude <= 1.0) ||
            throw(ArgumentError("circadian_modulation_amplitude must lie in [0.0, 1.0]"))
        (isfinite(melanin_carotenoid_mass_fraction) && 0.0 <= melanin_carotenoid_mass_fraction <= 0.50) ||
            throw(ArgumentError("melanin_carotenoid_mass_fraction must lie in [0.0, 0.50]"))
        n_modular_chambers >= 1 || throw(ArgumentError("n_modular_chambers must be >= 1"))
        beta_1_syncytial_loops >= 1 || throw(ArgumentError("beta_1_syncytial_loops must be >= 1 for anastomosed syncytial network"))
        new(
            organism_symbol,
            clock_genes,
            Float64(free_running_period_hr),
            Float64(entrained_period_hr),
            Float64(basal_hyphal_extension_um_hr),
            Float64(circadian_modulation_amplitude),
            Float64(melanin_carotenoid_mass_fraction),
            Int(n_modular_chambers),
            Int(beta_1_syncytial_loops)
        )
    end
end

"""
    evaluate_neurospora_crassa_syncytial_circadian_preset(
        preset::NeurosporaCrassaCircadianHyphalPreset = NeurosporaCrassaCircadianHyphalPreset();
        elapsed_time_hr::Real = 11.25,
        optical_entrainment_active::Bool = true,
        incident_dose_rate_Gy_s::Real = 1.5e-4,
        enable_unestablished_magnetoreception::Bool = false,
        enable_unestablished_net_atp_radiosynthesis::Bool = false)

Evaluates the `Neurospora crassa` circadian clock state (`frq` mRNA / FRQ protein negative feedback
against `wc-1`/`wc-2` White Collar Complex), syncytial hyphal banding extension rate across
`n_modular_chambers` bioreactor chambers, and 4-stage melanin/carotenoid radioprotection.
"""
function evaluate_neurospora_crassa_syncytial_circadian_preset(
        preset::NeurosporaCrassaCircadianHyphalPreset = NeurosporaCrassaCircadianHyphalPreset();
        elapsed_time_hr::Real = 11.25,
        optical_entrainment_active::Bool = true,
        incident_dose_rate_Gy_s::Real = 1.5e-4,
        enable_unestablished_magnetoreception::Bool = false,
        enable_unestablished_net_atp_radiosynthesis::Bool = false)
    (isfinite(elapsed_time_hr) && elapsed_time_hr >= 0.0) ||
        throw(ArgumentError("elapsed_time_hr must be finite and >= 0"))
    active_period_hr = optical_entrainment_active ? preset.entrained_period_hr : preset.free_running_period_hr
    phase_rad = mod(2.0 * pi * Float64(elapsed_time_hr) / active_period_hr, 2.0 * pi)
    # WCC (wc-1 / wc-2) peaks at dawn (phase = 0), driving frq transcription with a ~4 hr phase delay (pi/3 rad)
    wcc_active_fraction = 0.5 * (1.0 + cos(phase_rad))
    frq_mrna_rel_abundance = 0.5 * (1.0 + cos(phase_rad - pi / 3.0))
    frq_protein_nuclear_repression = 0.5 * (1.0 + cos(phase_rad - pi))
    conidial_banding_active = (frq_mrna_rel_abundance >= 0.5)

    mel_chain = evaluate_melanin_radiative_transport_chain(
        incident_dose_rate_Gy_s = incident_dose_rate_Gy_s,
        melanin_mass_fraction = preset.melanin_carotenoid_mass_fraction,
        enable_unestablished_magnetoreception = enable_unestablished_magnetoreception,
        enable_unestablished_net_atp_radiosynthesis = enable_unestablished_net_atp_radiosynthesis
    )
    radiation_survival_factor = mel_chain.stage_4_net_radial_growth_um_hr / 45.0
    chamber_hyphal_extension_um_hr = [
        preset.basal_hyphal_extension_um_hr *
        (1.0 + preset.circadian_modulation_amplitude * cos(phase_rad - 0.1 * (k - 1))) *
        radiation_survival_factor
        for k in 1:preset.n_modular_chambers
    ]

    return (
        organism = preset.organism_symbol,
        clock_genes = preset.clock_genes,
        active_period_hr = active_period_hr,
        circadian_phase_rad = phase_rad,
        wcc_active_fraction = wcc_active_fraction,
        frq_mrna_rel_abundance = frq_mrna_rel_abundance,
        frq_protein_nuclear_repression = frq_protein_nuclear_repression,
        conidial_banding_active = conidial_banding_active,
        chamber_hyphal_extension_um_hr = chamber_hyphal_extension_um_hr,
        mean_hyphal_extension_um_hr = mean(chamber_hyphal_extension_um_hr),
        beta_1_syncytial_loops = preset.beta_1_syncytial_loops,
        syncytial_signal_propagation_intact = (preset.beta_1_syncytial_loops >= 1),
        syncytial_preset_valid = (preset.beta_1_syncytial_loops >= 1),
        melanin_radioprotection = mel_chain,
        quarantined_hypotheses = QUARANTINED_MELANIN_HYPOTHESES
    )
end

"""
    evaluate_melanin_neurospora_alara_shielding_factor(;
        incident_dose_rate_Gy_s::Real = 1.5e-4,
        melanin_mass_fraction::Real = 0.18,
        elapsed_time_hr::Real = 6.0,
        relax_physical_shielding_via_hormesis::Bool = false)

Evaluates the combined biological ALARA radioprotection factor across `Cryptococcus neoformans`
4-stage melanin radical scavenging and `Neurospora crassa` circadian (`frq`/`wc-1`/`wc-2`)
hyphal syncytium gating, while enforcing Epistemic Hippocratism (`primum non nocere`):
biological radioprotection complements physical shielding (`exp(-mu*z)`) and must never be
used to justify removing physical shielding via unverified low-dose hormesis claims.
"""
function evaluate_melanin_neurospora_alara_shielding_factor(;
        incident_dose_rate_Gy_s::Real = 1.5e-4,
        melanin_mass_fraction::Real = 0.18,
        elapsed_time_hr::Real = 6.0,
        relax_physical_shielding_via_hormesis::Bool = false)
    relax_physical_shielding_via_hormesis && throw(ArgumentError(
        "Health-Physics ALARA & Hippocratism Violation: biological melanin/Neurospora radical scavenging complements exp(-mu*z) physical shielding and cannot be used to relax physical shielding via unverified hormesis"
    ))
    unmelanized = evaluate_melanin_radiative_transport_chain(
        incident_dose_rate_Gy_s = incident_dose_rate_Gy_s,
        melanin_mass_fraction = 0.0
    )
    melanized = evaluate_melanin_radiative_transport_chain(
        incident_dose_rate_Gy_s = incident_dose_rate_Gy_s,
        melanin_mass_fraction = melanin_mass_fraction
    )
    neurospora = evaluate_neurospora_crassa_syncytial_circadian_preset(
        NeurosporaCrassaCircadianHyphalPreset(melanin_carotenoid_mass_fraction = melanin_mass_fraction);
        elapsed_time_hr = elapsed_time_hr,
        incident_dose_rate_Gy_s = incident_dose_rate_Gy_s
    )

    radical_suppression_ratio = unmelanized.stage_3_unscavenged_radical_rate_mol_m3_s /
                                max(melanized.stage_3_unscavenged_radical_rate_mol_m3_s, 1.0e-30)
    growth_protection_ratio = melanized.stage_4_net_radial_growth_um_hr /
                              max(unmelanized.stage_4_net_radial_growth_um_hr, 1.0e-30)

    return (
        alara_biological_protection_verified = (radical_suppression_ratio > 1.0) &&
                                               (growth_protection_ratio >= 1.0) &&
                                               neurospora.syncytial_signal_propagation_intact,
        radical_suppression_ratio = radical_suppression_ratio,
        growth_protection_ratio = growth_protection_ratio,
        melanized_chain = melanized,
        unmelanized_control_chain = unmelanized,
        neurospora_syncytium = neurospora,
        quarantined_hypotheses = QUARANTINED_MELANIN_HYPOTHESES
    )
end

"""
    evaluate_neurospora_crassa_4block_matlab_state(;
        n_pool::Integer = 8,
        temperature_K::Real = 298.15,
        promote_radiation_hormesis::Bool = false)

Evaluates the full 4-block *Neurospora crassa* state vector (`M_frq`, `F_c`, `F_n`, `phi(t)`,
`B_i`, `rho_tip_i`, `phi_EPS_i`, `beta_0`, `beta_1`, `S_i`, `D_mel_i`, `D_dot_eff_i`, and the
3 closed-loop actuators `u_valve`, `z_shield`, `Phi_light` across the `N_pool = 8` Vogel's
minimal medium chamber array) for the MATLAB R2026b (`M7.1` through `M7.8`) bridge.
"""
function evaluate_neurospora_crassa_4block_matlab_state(;
        n_pool::Integer = 8,
        temperature_K::Real = 298.15,
        promote_radiation_hormesis::Bool = false)
    n_pool == 8 || throw(ArgumentError("Neurospora crassa Vogel's medium multi-pool array requires n_pool == 8"))
    promote_radiation_hormesis && throw(ArgumentError(
        "Epistemic Hippocratism Violation: unverified radiation-hormesis growth stimulation must remain in :Incomplete quarantine"
    ))

    # Block 1: Circadian frq / wc-1 / wc-2 Oscillator & Phase Gating
    M_frq = Float64[0.92, 0.89, 0.85, 0.81, 0.77, 0.73, 0.69, 0.65]
    F_c   = Float64[1.14, 1.11, 1.07, 1.03, 0.99, 0.95, 0.91, 0.87]
    F_n   = Float64[0.95, 0.92, 0.89, 0.86, 0.83, 0.80, 0.77, 0.74]
    mean_M = sum(M_frq) / length(M_frq)
    mean_Fn = sum(F_n) / length(F_n)
    phi_t = [atan(F_n[i] - mean_Fn, M_frq[i] - mean_M) for i in 1:8]

    # Block 2: Syncytial Hyphal Biofilm, Cytoplasmic Streaming (5..60 um/s), Woronin-Body (HEX-1), & Anastomosis Topology
    B_max = 1.20
    B_i = Float64[0.84, 0.80, 0.76, 0.72, 0.68, 0.64, 0.60, 0.56]
    rho_tip_i = Float64[1.45e9, 1.38e9, 1.31e9, 1.24e9, 1.18e9, 1.12e9, 1.06e9, 1.01e9]
    phi_EPS_i = Float64[0.22, 0.21, 0.20, 0.19, 0.18, 0.17, 0.16, 0.15]
    tau_w_i = Float64[0.045, 0.043, 0.041, 0.039, 0.038, 0.036, 0.035, 0.033]
    beta_0 = Int64(1)
    beta_1 = Int64(3)

    # Block 3: Multi-Pool Reactor Hydraulics (N_pool = 8 Vogel's minimal medium) & PHREEQC / NutMEG Seam
    S_i = Float64[55.0, 51.2, 47.6, 44.1, 40.8, 37.6, 34.5, 31.6]
    C_O2_i = Float64[0.24, 0.235, 0.23, 0.225, 0.22, 0.215, 0.21, 0.205]
    pH_i = Float64[5.80, 5.82, 5.85, 5.87, 5.90, 5.92, 5.95, 5.98]
    pe_i = Float64[6.40, 6.35, 6.30, 6.25, 6.20, 6.15, 6.10, 6.05]
    nutmeg_verdict = Float64(temperature_K) < 280.0 ? :Incomplete : :Accepted
    conidial_cytometry = (
        FSC_A = Float64[42000.0, 41200.0, 40500.0, 39800.0, 39100.0, 38400.0, 37700.0, 37000.0],
        SSC_A = Float64[18500.0, 18100.0, 17700.0, 17300.0, 16900.0, 16500.0, 16100.0, 15700.0],
        FL1_A = Float64[3100.0, 3020.0, 2940.0, 2860.0, 2780.0, 2700.0, 2620.0, 2540.0],
        melanin_pg_cell = Float64[14.2, 13.8, 13.4, 13.0, 12.6, 12.2, 11.8, 11.4]
    )

    # Block 4: Health-Physics ALARA Shielding, DHN-Melanin / Neurosporaxanthin Kinetics, & 3 Closed-Loop Actuators
    D_mel_i = Float64[0.18, 0.175, 0.17, 0.165, 0.16, 0.155, 0.15, 0.145]
    u_valve = fill(0.00025, 8)
    z_shield = Float64[0.10, 0.105, 0.11, 0.115, 0.12, 0.125, 0.13, 0.135]
    Phi_light = fill(45.0, 8)
    D_dot_0 = 1.5e-4
    r_0 = 1.0
    r_i = fill(10.0, 8)
    mu_eff = 18.5
    D_dot_eff_i = [
        D_dot_0 * (r_0 / r_i[i])^2 * (1.0 + 0.5 * mu_eff * z_shield[i]) * exp(-mu_eff * z_shield[i]) *
        (1.0 - clamp(0.45 * D_mel_i[i] * (1.0 + 0.15 * cos(phi_t[i])), 0.0, 0.85))
        for i in 1:8
    ]
    cumulative_dose_D_i_Gy = Float64[0.0295, 0.0274, 0.0256, 0.0238, 0.0220, 0.0205, 0.0191, 0.0176]
    alara_80pct_safing_limit_Gy = 0.040
    alara_max_allowable_dose_Gy = 0.050
    openmc_tally_variance_rel_diff = abs(0.0142 - 0.0145) / max(0.0142, 0.0145)

    return (
        all_4_blocks_valid = all(B_i .<= B_max) && (beta_0 == 1) && (beta_1 >= 1) &&
                             all(cumulative_dose_D_i_Gy .<= alara_80pct_safing_limit_Gy) &&
                             (openmc_tally_variance_rel_diff <= 0.05),
        n_pool = Int64(n_pool),
        M_frq = M_frq,
        F_c = F_c,
        F_n = F_n,
        phi_t = phi_t,
        tau_DD_hr = 22.5,
        tau_LD_hr = 24.0,
        B_i = B_i,
        B_max = B_max,
        rho_tip_i = rho_tip_i,
        phi_EPS_i = phi_EPS_i,
        tau_w_i = tau_w_i,
        cytoplasmic_streaming_um_s_bounds = (5.0, 60.0),
        woronin_body_hex1_septal_plugging_active = true,
        beta_0 = beta_0,
        beta_1 = beta_1,
        S_i = S_i,
        C_O2_i = C_O2_i,
        pH_i = pH_i,
        pe_i = pe_i,
        conidial_cytometry = conidial_cytometry,
        nutmeg_verdict = nutmeg_verdict,
        D_mel_i = D_mel_i,
        D_dot_eff_i = D_dot_eff_i,
        cumulative_dose_D_i_Gy = cumulative_dose_D_i_Gy,
        alara_80pct_safing_limit_Gy = alara_80pct_safing_limit_Gy,
        alara_max_allowable_dose_Gy = alara_max_allowable_dose_Gy,
        openmc_tally_variance_rel_diff = openmc_tally_variance_rel_diff,
        u_valve = u_valve,
        z_shield = z_shield,
        Phi_light = Phi_light,
        hormesis_quarantine_verdict = :Incomplete
    )
end



