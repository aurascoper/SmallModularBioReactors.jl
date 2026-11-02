# ==============================================================================
# src/contracts/SpatialFrame.jl
# Versioned Spatial Coordinate Frame Contract, Topology Gate, Reversible
# Payload-to-Reactor Coordinate Transform, and SHA-256 Frame Provenance
# ==============================================================================

const VALID_SMBR_MESH_KINDS = (
    :watertight_hexahedral_mesh,
    :conservative_voxel_grid,
    :gaussian_splat,
    :unstructured_point_cloud,
)

struct SMBRCoordinateTransform
    source_frame::Symbol
    target_frame::Symbol
    scale_to_meters::Float64
    rotation_matrix::Matrix{Float64}
    translation_m::Vector{Float64}

    function SMBRCoordinateTransform(;
            source_frame::Symbol = :reactor_local_um,
            target_frame::Symbol = :spacecraft_payload_m,
            scale_to_meters::Real = 1.0e-6,
            rotation_matrix::AbstractMatrix{<:Real} = Matrix{Float64}(I, 3, 3),
            translation_m::AbstractVector{<:Real} = zeros(Float64, 3))
        scale_to_meters > 0.0 || throw(ArgumentError("scale_to_meters must be strictly positive"))
        size(rotation_matrix) == (3, 3) || throw(ArgumentError("rotation_matrix must be 3x3"))
        length(translation_m) == 3 || throw(ArgumentError("translation_m must be a 3-vector"))
        R = Float64.(rotation_matrix)
        isapprox(R' * R, Matrix{Float64}(I, 3, 3); atol = 1.0e-9) ||
            throw(ArgumentError("rotation_matrix must be orthogonal (R' * R == I)"))
        isapprox(det(R), 1.0; atol = 1.0e-9) ||
            throw(ArgumentError("rotation_matrix must be right-handed proper rotation (det(R) == +1)"))
        new(source_frame, target_frame, Float64(scale_to_meters), R, Float64.(translation_m))
    end
end

function transform_point_forward(tf::SMBRCoordinateTransform, pt_source::AbstractVector{<:Real})
    length(pt_source) == 3 || throw(ArgumentError("pt_source must have length 3"))
    return tf.rotation_matrix * (Float64.(pt_source) .* tf.scale_to_meters) .+ tf.translation_m
end

function transform_point_inverse(tf::SMBRCoordinateTransform, pt_target_m::AbstractVector{<:Real})
    length(pt_target_m) == 3 || throw(ArgumentError("pt_target_m must have length 3"))
    return (tf.rotation_matrix' * (Float64.(pt_target_m) .- tf.translation_m)) ./ tf.scale_to_meters
end

struct SMBRSpatialFrame
    schema_version::String
    frame_id::String
    mesh_kind::Symbol
    is_watertight_manifold::Bool
    closed_face_fluxes::Bool
    cell_count::Int
    cell_volumes_m3::Vector{Float64}
    cell_centroids_m::Matrix{Float64}
    beta_0::Int
    beta_1::Int
    beta_2::Int
    frame_sha256::String

    function SMBRSpatialFrame(;
            schema_version::AbstractString = "1.0.0",
            frame_id::AbstractString = "smbr_default_4cell",
            mesh_kind::Symbol = :watertight_hexahedral_mesh,
            is_watertight_manifold::Bool = true,
            closed_face_fluxes::Bool = true,
            cell_volumes_m3::AbstractVector{<:Real} = fill(1.0e-9, 4),
            cell_centroids_m::AbstractMatrix{<:Real} = [0.5e-3 1.5e-3 2.5e-3 3.5e-3; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0],
            beta_0::Int = 1,
            beta_1::Int = 2,
            beta_2::Int = 0)
        mesh_kind in VALID_SMBR_MESH_KINDS ||
            throw(ArgumentError("Unsupported mesh_kind :$(mesh_kind)"))
        n_cells = length(cell_volumes_m3)
        n_cells >= 1 || throw(ArgumentError("cell_volumes_m3 must contain at least 1 cell"))
        all(v -> isfinite(v) && v > 0.0, cell_volumes_m3) ||
            throw(ArgumentError("All cell_volumes_m3 must be finite and strictly positive"))
        size(cell_centroids_m) == (3, n_cells) ||
            throw(ArgumentError("cell_centroids_m must have size (3, $(n_cells))"))
        (beta_0 >= 0 && beta_1 >= 0 && beta_2 >= 0) ||
            throw(ArgumentError("Betti numbers must be non-negative"))

        vols = Float64.(cell_volumes_m3)
        cents = Float64.(cell_centroids_m)
        io = IOBuffer()
        print(io, schema_version, "|", frame_id, "|", mesh_kind, "|", is_watertight_manifold, "|", closed_face_fluxes, "|")
        for v in vols
            @printf(io, "%.12e,", v)
        end
        for c in cents
            @printf(io, "%.12e,", c)
        end
        print(io, "|", beta_0, ":", beta_1, ":", beta_2)
        digest = "sha256:" * bytes2hex(sha256(take!(io)))

        new(
            String(schema_version),
            String(frame_id),
            mesh_kind,
            is_watertight_manifold,
            closed_face_fluxes,
            n_cells,
            vols,
            cents,
            beta_0,
            beta_1,
            beta_2,
            digest
        )
    end
end

"""
    validate_spatial_frame_for_fv_solver(frame::SMBRSpatialFrame; require_conservative_mesh::Bool = true)

Verifies that a spatial frame has an intact SHA-256 digest and, when `require_conservative_mesh == true`,
rejects non-conservative point clouds or Gaussian splats from participating in finite-volume PDE flux solves.
"""
function validate_spatial_frame_for_fv_solver(frame::SMBRSpatialFrame; require_conservative_mesh::Bool = true)
    rebuilt = SMBRSpatialFrame(
        schema_version = frame.schema_version,
        frame_id = frame.frame_id,
        mesh_kind = frame.mesh_kind,
        is_watertight_manifold = frame.is_watertight_manifold,
        closed_face_fluxes = frame.closed_face_fluxes,
        cell_volumes_m3 = frame.cell_volumes_m3,
        cell_centroids_m = frame.cell_centroids_m,
        beta_0 = frame.beta_0,
        beta_1 = frame.beta_1,
        beta_2 = frame.beta_2
    )
    rebuilt.frame_sha256 == frame.frame_sha256 ||
        throw(ArgumentError("SMBRSpatialFrame SHA-256 integrity check failed for $(frame.frame_id)"))

    fv_eligible = (frame.mesh_kind in (:watertight_hexahedral_mesh, :conservative_voxel_grid)) &&
                  frame.is_watertight_manifold &&
                  frame.closed_face_fluxes
    if require_conservative_mesh && !fv_eligible
        throw(ArgumentError(
            "Spatial frame $(frame.frame_id) with mesh_kind=:$(frame.mesh_kind) is not eligible for conservative finite-volume PDE solves"
        ))
    end
    return (
        valid = true,
        fv_eligible = fv_eligible,
        total_volume_m3 = sum(frame.cell_volumes_m3),
        euler_characteristic = frame.beta_0 - frame.beta_1 + frame.beta_2,
        frame_sha256 = frame.frame_sha256
    )
end
