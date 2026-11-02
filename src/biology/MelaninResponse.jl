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
