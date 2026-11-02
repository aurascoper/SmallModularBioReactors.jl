# ==============================================================================
# src/transport/PolymerTransport.jl
# Extracellular Polymeric Substance (EPS) & Biopolymer Transport, Cross-Linking,
# and Microfluidic Micropillar Streamer Time-Series (Demonstration 1: t = 5..13 hr)
# ==============================================================================

const MICROPILLAR_STREAMER_TIMEPOINTS_HR = (5.0, 7.0, 9.0, 11.0, 12.0, 13.0)

"""
    simulate_micropillar_streamer_bridging(;
        timepoints_hr::AbstractVector{<:Real} = collect(MICROPILLAR_STREAMER_TIMEPOINTS_HR),
        initial_eps_kg_m3::Real = 0.08,
        initial_porosity::Real = 0.42,
        shear_rate_s_inv::Real = 15.0,
        divalent_ca_mol_m3::Real = 1.5,
        pillar_pitch_um::Real = 20.0)

Simulates green-fluorescent biofilm EPS streamer elongation, divalent-cation cross-linking,
and pore-bridging anastomosis across circular micropillars (`20 um` scale bar, `t = 5, 7, 9, 11, 12, 13 hr`
from Demonstration 1 image `rOkmBjSK...jpg`).
"""
function simulate_micropillar_streamer_bridging(;
        timepoints_hr::AbstractVector{<:Real} = collect(MICROPILLAR_STREAMER_TIMEPOINTS_HR),
        initial_eps_kg_m3::Real = 0.08,
        initial_porosity::Real = 0.42,
        shear_rate_s_inv::Real = 15.0,
        divalent_ca_mol_m3::Real = 1.5,
        pillar_pitch_um::Real = 20.0)
    isempty(timepoints_hr) && throw(ArgumentError("timepoints_hr cannot be empty"))
    issorted(timepoints_hr) || throw(ArgumentError("timepoints_hr must be sorted ascending"))
    (initial_eps_kg_m3 > 0.0 && 0.05 < initial_porosity < 0.95 && shear_rate_s_inv > 0.0 &&
     divalent_ca_mol_m3 >= 0.0 && pillar_pitch_um > 0.0) ||
        throw(ArgumentError("Invalid physical parameters for micropillar streamer bridging"))

    n_t = length(timepoints_hr)
    eps_trajectory = zeros(Float64, n_t)
    crosslink_fraction = zeros(Float64, n_t)
    porosity_trajectory = zeros(Float64, n_t)
    permeability_m2 = zeros(Float64, n_t)
    beta_1_loops = zeros(Int, n_t)

    k0 = 3.5e-13
    for idx in 1:n_t
        t_hr = Float64(timepoints_hr[idx])
        dt_from_5 = max(0.0, t_hr - 5.0)
        crosslink_fraction[idx] = clamp(
            (Float64(divalent_ca_mol_m3) / (0.5 + Float64(divalent_ca_mol_m3))) * (1.0 - exp(-0.22 * t_hr)),
            0.0,
            0.98
        )
        eps_val = Float64(initial_eps_kg_m3) * (1.0 + 0.32 * dt_from_5 * (1.0 + 0.4 * crosslink_fraction[idx]))
        eps_trajectory[idx] = eps_val
        clog_frac = clamp(0.18 * (eps_val - Float64(initial_eps_kg_m3)), 0.0, Float64(initial_porosity) - 0.08)
        phi_t = Float64(initial_porosity) - clog_frac
        porosity_trajectory[idx] = phi_t
        permeability_m2[idx] = compute_kozeny_carman_permeability(phi_t; phi0 = initial_porosity, k0_m2 = k0)
        # Anastomosis pore-bridging loop closure across micropillars as t progresses from 5 to 13 hr
        beta_1_loops[idx] = t_hr < 7.0 ? 0 : (t_hr < 11.0 ? 1 : 2)
    end

    return (
        timepoints_hr = Float64.(timepoints_hr),
        pillar_scale_bar_um = Float64(pillar_pitch_um),
        eps_kg_m3 = eps_trajectory,
        crosslink_fraction = crosslink_fraction,
        porosity = porosity_trajectory,
        permeability_m2 = permeability_m2,
        beta_1_anastomosis_loops = beta_1_loops,
        streamer_bridging_established = (last(beta_1_loops) >= 2 && last(permeability_m2) < first(permeability_m2)),
        evidence_image_sha256 = DOWNLOADS_FIVE_IMAGES_MANIFEST.demonstration_1_microfluidic_biofilm_streamers.sha256
    )
end
