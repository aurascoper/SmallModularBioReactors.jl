# ==============================================================================
# src/molecular/SobolevAdapter.jl
# SobolevMacromolecule (SMACRO1 Binary Payload) & Molecular Anchor Adapter
# Strictly Separating Discrete Sobolev Relaxation Steps from Physical Seconds
# ==============================================================================

const SMACRO1_MAGIC_BYTES = UInt8['S', 'M', 'A', 'C', 'R', 'O', '1', 0x00]

struct SobolevMolecularAnchor
    anchor_id::String
    cell_index::Int
    sobolev_relaxation_step::Int
    bead_coordinates_nm::Matrix{Float64}
    payload_sha256::String

    function SobolevMolecularAnchor(;
            anchor_id::AbstractString = "csgA_curli_anchor_1",
            cell_index::Int = 1,
            sobolev_relaxation_step::Int = 100,
            bead_coordinates_nm::AbstractMatrix{<:Real} = [0.0 0.38 0.76 1.14; 0.0 0.12 -0.10 0.08; 0.0 0.05 0.10 0.15])
        isempty(strip(anchor_id)) && throw(ArgumentError("anchor_id cannot be empty"))
        cell_index >= 1 || throw(ArgumentError("cell_index must be >= 1"))
        sobolev_relaxation_step >= 0 || throw(ArgumentError("sobolev_relaxation_step must be >= 0"))
        (size(bead_coordinates_nm, 1) == 3 && size(bead_coordinates_nm, 2) >= 1) ||
            throw(ArgumentError("bead_coordinates_nm must have size (3, N_beads >= 1)"))
        all(isfinite, bead_coordinates_nm) || throw(ArgumentError("bead_coordinates_nm must be finite"))

        coords = Float64.(bead_coordinates_nm)
        raw = encode_smbr_smacro1_bytes(sobolev_relaxation_step, coords)
        digest = "sha256:" * bytes2hex(sha256(raw))
        new(String(anchor_id), cell_index, sobolev_relaxation_step, coords, digest)
    end
end

"""
    encode_smbr_smacro1_bytes(sobolev_step::Int, bead_coords_nm::AbstractMatrix{<:Real}) -> Vector{UInt8}

Encodes a binary `SMACRO1` frame (8-byte magic header + Int32 step + Int32 bead count + Float64 xyz triples).
"""
function encode_smbr_smacro1_bytes(sobolev_step::Int, bead_coords_nm::AbstractMatrix{<:Real})
    sobolev_step >= 0 || throw(ArgumentError("sobolev_step must be >= 0"))
    size(bead_coords_nm, 1) == 3 || throw(ArgumentError("bead_coords_nm must be 3xN"))
    n_beads = size(bead_coords_nm, 2)
    io = IOBuffer()
    write(io, SMACRO1_MAGIC_BYTES)
    write(io, Int32(sobolev_step))
    write(io, Int32(n_beads))
    for j in 1:n_beads, i in 1:3
        write(io, Float64(bead_coords_nm[i, j]))
    end
    return take!(io)
end

"""
    decode_smbr_smacro1_bytes(bytes::AbstractVector{UInt8})

Decodes and validates a binary `SMACRO1` frame, throwing `ArgumentError` on corrupted headers
or truncated coordinate payloads.
"""
function decode_smbr_smacro1_bytes(bytes::AbstractVector{UInt8})
    length(bytes) >= 16 || throw(ArgumentError("SMACRO1 payload truncated (< 16 bytes)"))
    bytes[1:8] == SMACRO1_MAGIC_BYTES || throw(ArgumentError("Invalid SMACRO1 magic header"))
    io = IOBuffer(Vector{UInt8}(bytes))
    read(io, 8)
    step = Int(read(io, Int32))
    n_beads = Int(read(io, Int32))
    (step >= 0 && n_beads >= 1) || throw(ArgumentError("Invalid SMACRO1 step or bead count"))
    expected_len = 16 + 24 * n_beads
    length(bytes) == expected_len || throw(ArgumentError("SMACRO1 byte length mismatch: expected $(expected_len), got $(length(bytes))"))
    coords = zeros(Float64, 3, n_beads)
    for j in 1:n_beads, i in 1:3
        coords[i, j] = read(io, Float64)
    end
    return (
        sobolev_relaxation_step = step,
        bead_count = n_beads,
        bead_coordinates_nm = coords
    )
end
