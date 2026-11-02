# ==============================================================================
# src/inference/EvidenceGates.jl
# Forwarding bridge to src/verification/EvidenceGates.jl so callers referencing
# src/inference/EvidenceGates.jl resolve the canonical NeuroARC 3-Valued and
# Epistemic Hippocratism (Primum Non Nocere) quarantine gates without duplication.
# Strictly zero em-dashes and zero en-dashes.
# ==============================================================================

if !isdefined(@__MODULE__, :evaluate_neuroarc_evidence_gate)
    include(joinpath(@__DIR__, "..", "verification", "EvidenceGates.jl"))
end
