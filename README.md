# SmallModularBioReactors.jl

`SmallModularBioReactors.jl` is a bounded, read-only, test-first Julia orchestration and verification package for modular interplanetary biological and biochemical reactors (`NeuralCompose x M7 x SPICE.jl x Basilisk`).

It enforces the strict architectural separation:
`Physical solver != Inference engine != Evidence verifier != Visualization`
and explicitly distinguishes **Small Modular Biological Reactors** (e.g. Sartorius `Ambr 250 Modular` benchtop bioprocess stations and ESA `BIM` microfluidic cassettes) from external **Small Modular Fission Power Systems** (e.g. Mars `Kilopower` surface electrical/thermal power busses). It is strictly a biological/biochemical bioprocess engineering and reactive-transport verification layer, NOT a nuclear fission criticality simulator.

## 1. Three-System Architecture

1. **Biological and Geochemical Reactor (`src/contracts/`, `src/transport/`, `src/biology/`, `src/molecular/`)**:
   - Coupled state vector $X = [u, p, T, \phi, k, C_{\mathrm{aq}}, M_{\mathrm{min}}, N_{\mathrm{iso}}, B, E, P]$ with global elemental conservation $dM_i^{\mathrm{total}}/dt = J_{i,\mathrm{in}} - J_{i,\mathrm{out}} + R_i^{\mathrm{nuclear}}$.
   - 3-fraction analytical accounting (`:ICP_MS`, `:XRD_XRF`, `:Flow_Cytometry`) and System A (`Sudoite + illite`), System B (`Bentonite / montmorillonite`), and System C (`Illite/smectite + EPS/biofilm`) clay transitions.
   - `Menke et al. (2026)` mass-conservative Darcy-Brinkman-Stokes (DBS) reactive transport (`GeoChemFoam` + `PHREEQC`) coupling mobile pore-water advection/diffusion with cell-scale adsorbed surface inventories (`P-L + M^{2+} <=> P-L-M`, mass residual $\ll 1\%$).
   - Conservative micro-continuum-cellular automaton (CA) biofilm redistribution (`arXiv:2610.01565`), `Neurospora crassa` syncytial circadian (`frq`/`wc-1`/`wc-2`) and melanin/carotenoid radioprotection preset (`NeurosporaCrassaCircadianHyphalPreset`, `evaluate_neurospora_crassa_syncytial_circadian_preset`), `:Arc` capsid (`PDB:6GSE`) and peptide/protein registry (`CsgA`, `DsrB`, `MtrC`, `HydA`, `LanM`), and the **Section 5 Traceable 7 Budget-Closure Criteria** (`NutMEG` Higgins & Cockell 2020 PMC7653372 with mandatory $<280\text{ K}$ out-of-range flag, `Affholder et al. 2021` Enceladus plume Bayes factor, `Sauterey et al. 2022` early Mars, `Europa` held in `:Incomplete`, and `PlanetProfile` Styczinski et al. 2023 SPICE-to-interior coupling).

2. **Spacecraft and Planetary-Environment Digital Twin (`src/environment/`, `src/spacecraft/`, `src/visualization/`)**:
   - `SPICEAdapter.jl`: `SPICE.jl` (`v0.2.5`, Julia `>= 1.6`, `CSPICE`) and `Basilisk` `spiceInterface` (`SpicePlanetStateMsgPayload` in SI `m` and `m/s` with `J20002Pfix` `IAU_` frames) connected via the single boundary converter `convert_spice_state_to_basilisk_si(pos_km, vel_km_s, dcm_j2000_to_iau)`, serialized under `SMBR_SPICE_LOCK = ReentrantLock()`, and strictly separated from `RadiationField.jl`.
   - `VizardAdapter.jl`: Unity-based `Vizard` (`vizInterface` protobuf over ZeroMQ with handshake strings `PING`, `SIM_UPDATE`, `OK`, `REQUEST_INPUT`, `VIZARD_INPUT`, `REQUEST_IMAGE_`) supporting `liveStream` (`tcp://localhost:5556`), `broadcastStream` (`tcp://localhost:5570` non-blocking read-only HUD tap), and `noDisplay`, paired with `M7Adapter.jl` across the Three-Clock Law (`MissionClock.jl`).

3. **Scientific Inference and Evidence Verification (`src/inference/`, `src/verification/`)**:
   - `BayesianCalibration.jl` ($p(\theta \mid y) \propto p(y \mid \theta)p(\theta)$) and `HastingsMetropolis.jl` ($\alpha(\theta, \theta')$, ESS, acceptance rate, autocorrelation, PPC, strictly separated from engineering build decisions).
   - `DriftDetection.jl` ($\Delta_t = \mu_{\theta,t} - \mu_{\theta,t-1}$ and $\Delta\Sigma_t = \Sigma_{\theta,t} - \Sigma_{\theta,t-1}$ stored separately from physical state change, priors, and observation models) and 6-objective `ParetoDesign.jl` (`arXiv:2508.10970`).
   - `EvidenceGates.jl`: NeuroARC 3-valued verdicts (`:Accepted`, `:Rejected`, `:Incomplete`) where `:Incomplete` never collapses into `:Rejected` without explicit refutation evidence.

## 2. Five `~/Downloads` Visual Evidence Anchors (`2026-11-01 23:34..23:35`)

- **Demonstration 1 (`rOkmBjSK...jpg`, `899x553`, `75,741 bytes`)**: Fluorescence microscopy time-series (`t = 5, 7, 9, 11, 12, 13 hr`, `20 um` scale bar) of green-fluorescent biofilm streamer formation and pore-bridging anastomosis across circular micropillars.
- **Demonstration 2 (`4I3OYyXt...jpg`, `736x652`, `20,116 bytes`)**: Paired fungal colony phenotypes on agar plate comparing heavily melanized dark radial hyphae (left) vs lighter radially symmetric control colony (right).
- **Demonstration 3 (`ytJAGCW0...jpg`, `960x834`, `137,863 bytes`)**: CAD assembly of the ESA Biology in Microgravity (`BIM`) cylindrical carousel payload module with modular sample cassettes and central rotor.
- **Section 9 External Power Boundary (`MzHXxnc_...jpg`, `1200x801`, `65,914 bytes`)**: Planetary surface lander and four Kilopower-style umbrella radiator fission power units tethered on Mars regolith (external electrical/thermal power bus only).
- **Section 9 Bioprocess Workstation (`aqLGjX4F...jpg`, `800x600`, `27,715 bytes`)**: Sartorius `Ambr 250 Modular` multi-parallel benchtop bioreactor workstation.

## 3. Running the Test Suite

```bash
julia --project=. test/runtests.jl
```
