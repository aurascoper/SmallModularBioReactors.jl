# ==============================================================================
# src/molecular/StructureProvenance.jl
# Strict Separation Between 3D Coordinate Inspection (PyMOL / UCSF ChimeraX)
# and Experimentally Measured Thermodynamic Binding Constants (K_d / K_bind)
# ==============================================================================

struct StructureAndBindingProvenance
    protein_symbol::Symbol
    coordinate_viewer::Symbol
    pdb_or_model_tag::String
    is_experimental_coordinates::Bool
    binding_constant_K_bind_m3_mol::Union{Nothing,Float64}
    binding_assay_method::Union{Nothing,Symbol}
    provenance_separated::Bool

    function StructureAndBindingProvenance(;
            protein_symbol::Symbol = :CsgA,
            coordinate_viewer::Symbol = :PyMOL,
            pdb_or_model_tag::AbstractString = "[THEORETICAL_MODEL_FRAGMENT]",
            is_experimental_coordinates::Bool = false,
            binding_constant_K_bind_m3_mol::Union{Nothing,Real} = 25.0,
            binding_assay_method::Union{Nothing,Symbol} = :ITC_and_ICP_MS_isotherm,
            claim_kd_from_static_pdb_coordinates_only::Bool = false)
        coordinate_viewer in (:PyMOL, :ChimeraX, :SobolevMacromolecule) ||
            throw(ArgumentError("Unsupported coordinate_viewer :$(coordinate_viewer)"))
        claim_kd_from_static_pdb_coordinates_only && throw(ArgumentError(
            "Provenance Separation Violation: Thermodynamic binding constants (K_d / K_bind) cannot be derived from static PyMOL/ChimeraX coordinate inspection alone; an experimental assay method (e.g. :ITC_and_ICP_MS_isotherm, :SPR) is required"
        ))
        if binding_constant_K_bind_m3_mol !== nothing
            binding_constant_K_bind_m3_mol > 0.0 || throw(ArgumentError("binding_constant_K_bind_m3_mol must be > 0"))
            binding_assay_method !== nothing || throw(ArgumentError(
                "Providing a numeric binding_constant_K_bind_m3_mol requires an explicit binding_assay_method"
            ))
        end

        new(
            protein_symbol,
            coordinate_viewer,
            String(pdb_or_model_tag),
            is_experimental_coordinates,
            binding_constant_K_bind_m3_mol === nothing ? nothing : Float64(binding_constant_K_bind_m3_mol),
            binding_assay_method,
            true
        )
    end
end
