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
        melanin_radioprotection = mel_chain,
        quarantined_hypotheses = QUARANTINED_MELANIN_HYPOTHESES
    )
end

