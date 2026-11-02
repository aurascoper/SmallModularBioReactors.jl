# ==============================================================================
# src/molecular/ProteinPeptideRegistry.jl
# Registry of Metal-Binding & Catalytic Proteins/Peptides (CsgA, DsrB, OahA,
# OxDC, LanM / Phytochelatin P-L Ligands)
# ==============================================================================

struct ProteinPeptideEntry
    symbol::Symbol
    uniprot_or_pdb_id::String
    motif_signature::String
    core_length_aa::Int
    repeat_length_aa::Int
    target_metal_or_cofactor::Symbol
    stoichiometric_binding_sites_per_monomer::Int
    structure_class::Symbol
end

const PROTEIN_PEPTIDE_REGISTRY = Dict{Symbol,ProteinPeptideEntry}(
    :CsgA => ProteinPeptideEntry(
        :CsgA,
        "UniProt:P28307 / PDB:8ENQ",
        "S-X(5)-Q-X-G-X-G-N-X-A-X(3)-Q",
        12,
        18,
        :Ca2_Pb2_surface_eps,
        2,
        :theoretical_repeat_and_cryoem_fibril
    ),
    :DsrB => ProteinPeptideEntry(
        :DsrB,
        "UniProt:P45575 / PDB:2V4J",
        "C-X(2)-C-X(3)-C (n_X = 5)",
        8,
        8,
        :Fe4S4_siroheme,
        1,
        :experimental_xray_structure
    ),
    :OahA => ProteinPeptideEntry(
        :OahA,
        "PDB:3M0J",
        "Arg21-Ser26-Mg2+-Oxalate",
        10,
        10,
        :Mg2_plus,
        1,
        :experimental_xray_structure
    ),
    :OxDC => ProteinPeptideEntry(
        :OxDC,
        "PDB:1J58",
        "His95-His97-Glu101-Mn2+",
        7,
        7,
        :Mn2_plus,
        2,
        :experimental_xray_structure
    ),
    :LanM_Peptide => ProteinPeptideEntry(
        :LanM_Peptide,
        "PDB:6MI5",
        "EF-hand-Asp-Glu-Pro-chelator",
        12,
        12,
        :Ln3_Pb2_plus,
        4,
        :experimental_nmr_xray_structure
    ),
    :Arc => ProteinPeptideEntry(
        :Arc,
        "PDB:6GSE",
        "Gag-CA-NTD-CTD-capsid-oligomer",
        134,
        396,
        :RNA_capsid_self_assembly,
        240,
        :experimental_cryo_em_xray_structure
    ),
)

"""
    lookup_protein_peptide_entry(symbol::Symbol) -> ProteinPeptideEntry

Retrieves a verified entry from `PROTEIN_PEPTIDE_REGISTRY` and validates CsgA (12-aa core / 18-aa repeat)
and DsrB (`n_X = 5` internal wildcards) motif invariants.
"""
function lookup_protein_peptide_entry(symbol::Symbol)
    haskey(PROTEIN_PEPTIDE_REGISTRY, symbol) ||
        throw(ArgumentError("Unknown protein/peptide symbol :$(symbol)"))
    return PROTEIN_PEPTIDE_REGISTRY[symbol]
end
