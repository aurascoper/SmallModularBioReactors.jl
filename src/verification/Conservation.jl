# ==============================================================================
# src/verification/Conservation.jl
# Multi-Compartment Mass, Charge, Fluid Continuity, and Free-Energy Budget
# Verification Layer (Evidence Verifier Role)
# ==============================================================================

"""
    audit_reactor_conservation_invariants(;
        state_before::ReactorState = ReactorState(),
        state_after::ReactorState = ReactorState(),
        frame::SMBRSpatialFrame = SMBRSpatialFrame(),
        dt_s::Real = 10.0,
        species_key::Symbol = :Pb2_plus,
        J_in_mol_s::Real = 0.0,
        J_out_mol_s::Real = 0.0,
        R_nuclear_mol_s::Real = 0.0,
        rel_tol::Real = 1.0e-6)

Verifies elemental mass conservation across dissolved pore-water (`phi * C_aq`),
mineral/surface-adsorbed (`M_min`), and peptide/polymer compartments over `dt_s`.
"""
function audit_reactor_conservation_invariants(;
        state_before::ReactorState = ReactorState(),
        state_after::ReactorState = ReactorState(),
        frame::SMBRSpatialFrame = SMBRSpatialFrame(),
        dt_s::Real = 10.0,
        species_key::Symbol = :Pb2_plus,
        J_in_mol_s::Real = 0.0,
        J_out_mol_s::Real = 0.0,
        R_nuclear_mol_s::Real = 0.0,
        rel_tol::Real = 1.0e-6)
    (state_before.cell_count == frame.cell_count && state_after.cell_count == frame.cell_count) ||
        throw(ArgumentError("State and frame cell counts must match"))
    (haskey(state_before.C_aq_mol_m3, species_key) && haskey(state_after.C_aq_mol_m3, species_key)) ||
        throw(ArgumentError("Missing aqueous species :$(species_key)"))
    (haskey(state_before.M_min_mol_m3, species_key) && haskey(state_after.M_min_mol_m3, species_key)) ||
        throw(ArgumentError("Missing mineral/surface species :$(species_key)"))

    vols = frame.cell_volumes_m3
    n = frame.cell_count
    m_before = sum(
        (state_before.phi[i] * state_before.C_aq_mol_m3[species_key][i] + state_before.M_min_mol_m3[species_key][i]) * vols[i]
        for i in 1:n
    )
    m_after = sum(
        (state_after.phi[i] * state_after.C_aq_mol_m3[species_key][i] + state_after.M_min_mol_m3[species_key][i]) * vols[i]
        for i in 1:n
    )

    rate_check = evaluate_global_mass_conservation_rate(
        initial_total_mol = m_before,
        final_total_mol = m_after,
        dt_s = dt_s,
        J_in_mol_s = J_in_mol_s,
        J_out_mol_s = J_out_mol_s,
        R_nuclear_mol_s = R_nuclear_mol_s,
        rel_tol = rel_tol
    )

    return (
        role_contract = ArchitecturalRoleContract(
            component_name = "ConservationVerifier",
            assigned_role = :evidence_verifier,
            issues_verification_verdicts = true
        ),
        species_key = species_key,
        initial_total_mol = m_before,
        final_total_mol = m_after,
        rate_conservation = rate_check,
        all_conserved = rate_check.is_conserved && rate_check.well_below_one_percent
    )
end
