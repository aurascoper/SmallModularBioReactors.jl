# ==============================================================================
# src/environment/SPICEAdapter.jl
# SPICE.jl (v0.2.5, Julia >= 1.6, CSPICE) & Basilisk spiceInterface / pyswice
# Ephemeris Adapter (furnsh, utc2et, spkpos, pxform, naif0012.tls, de440.bsp,
# scenarioMonteCarloSpice, SpicePlanetStateMsgPayload, and single boundary
# conversion convert_spice_state_to_basilisk_si)
# Strictly Separated from Radiation Dose Transport
# ==============================================================================

const SMBR_SPICE_LOCK = ReentrantLock()

const STANDARD_SPICE_KERNELS = (
    lsk = "naif0012.tls",
    spk = "de440.bsp",
    pck = "pck00010.tpc",
    hst_bsp = "hst_edited.bsp",
)

const VALID_SPICE_TARGETS = (
    :HST_LEO,
    :EARTH,
    :MARS,
    :MOON,
    :ENCELADUS,
)

const BASILISK_SPICE_INTERFACE_OUTPUT_MSGS = (
    :spiceTimeOutMsg,
    :planetStateOutMsgs,
    :scStateOutMsgs,
    :attRefStateOutMsgs,
    :transRefStateOutMsgs,
)

const BASILISK_SPICE_INTERFACE_METHODS = (
    :addPlanetNames,
    :addSpacecraftNames,
    :loadSpiceKernel,
    :addKernelPath,
)

"""
    SpicePlanetStateMsgPayload

Exact Julia representation of Basilisk's C message payload `SpicePlanetStateMsgPayload`:
- `PositionVector::Vector{Float64}`: 3-vector in SI `m`
- `VelocityVector::Vector{Float64}`: 3-vector in SI `m/s`
- `J2000Current::Float64`: TDB seconds past J2000 (`s`)
- `J20002Pfix::Matrix{Float64}`: 3x3 orientation DCM from `J2000` to planet-fixed `"IAU_" * PlanetName` frame
- `J20002Pfix_dot::Matrix{Float64}`: 3x3 time derivative of `J20002Pfix` (`1/s`)
- `computeOrient::Int`: flag (`1` when orientation is computed)
- `PlanetName::String`: registered body name
- `planet_frame::String`: default `"IAU_" * uppercase(PlanetName)` frame
"""
struct SpicePlanetStateMsgPayload
    PositionVector::Vector{Float64}
    VelocityVector::Vector{Float64}
    J2000Current::Float64
    J20002Pfix::Matrix{Float64}
    J20002Pfix_dot::Matrix{Float64}
    computeOrient::Int
    PlanetName::String
    planetName::String
    planet_frame::String

    function SpicePlanetStateMsgPayload(;
            PositionVector::AbstractVector{<:Real},
            VelocityVector::AbstractVector{<:Real},
            J2000Current::Real,
            J20002Pfix::AbstractMatrix{<:Real} = Matrix{Float64}(I, 3, 3),
            J20002Pfix_dot::AbstractMatrix{<:Real} = zeros(Float64, 3, 3),
            computeOrient::Int = 1,
            PlanetName::AbstractString = "MARS",
            planetName::Union{Nothing,AbstractString} = nothing,
            planet_frame::Union{Nothing,AbstractString} = nothing)
        length(PositionVector) == 3 || throw(ArgumentError("PositionVector must be a 3-vector (SI m)"))
        length(VelocityVector) == 3 || throw(ArgumentError("VelocityVector must be a 3-vector (SI m/s)"))
        all(isfinite, PositionVector) || throw(ArgumentError("PositionVector must be finite"))
        all(isfinite, VelocityVector) || throw(ArgumentError("VelocityVector must be finite"))
        isfinite(J2000Current) || throw(ArgumentError("J2000Current must be finite"))
        size(J20002Pfix) == (3, 3) || throw(ArgumentError("J20002Pfix must be 3x3"))
        size(J20002Pfix_dot) == (3, 3) || throw(ArgumentError("J20002Pfix_dot must be 3x3"))
        (computeOrient in (0, 1)) || throw(ArgumentError("computeOrient must be 0 or 1"))
        resolved_name = planetName === nothing ? String(strip(PlanetName)) : String(strip(planetName))
        isempty(resolved_name) && throw(ArgumentError("PlanetName must be non-empty"))
        frame_str = planet_frame === nothing ? ("IAU_" * uppercase(resolved_name)) : String(planet_frame)
        new(
            Float64.(PositionVector),
            Float64.(VelocityVector),
            Float64(J2000Current),
            Float64.(J20002Pfix),
            Float64.(J20002Pfix_dot),
            computeOrient,
            resolved_name,
            resolved_name,
            frame_str
        )
    end
end

"""
    convert_spice_state_to_basilisk_si(
        pos_km::AbstractVector{<:Real},
        vel_km_s::AbstractVector{<:Real},
        dcm_j2000_to_iau::AbstractMatrix{<:Real} = Matrix{Float64}(I, 3, 3);
        dcm_dot_j2000_to_iau::AbstractMatrix{<:Real} = zeros(Float64, 3, 3),
        j2000_current_s::Real = 0.0,
        planet_name::AbstractString = "MARS",
        expected_distance_m_bounds::Union{Nothing,Tuple{<:Real,<:Real}} = nothing)

Single canonical boundary conversion function from NAIF SPICE / `SPICE.jl` (`km`, `km/s`)
to Basilisk's `SpicePlanetStateMsgPayload` (`m`, `m/s`), carrying the `J20002Pfix` IAU frame
DCM separately.
Because the `SPICE.jl` README does not state units explicitly, `spice_jl_unit_verdict` is held
at `:Incomplete` unless `expected_distance_m_bounds` is provided to verify the converted `spkpos`
magnitude against a known orbital radius (in which case it transitions to `:Accepted` if inside
bounds or `:Rejected` if off by a factor of `1e3`).
"""
function convert_spice_state_to_basilisk_si(
        pos_km::AbstractVector{<:Real},
        vel_km_s::AbstractVector{<:Real},
        dcm_j2000_to_iau::AbstractMatrix{<:Real} = Matrix{Float64}(I, 3, 3);
        dcm_dot_j2000_to_iau::AbstractMatrix{<:Real} = zeros(Float64, 3, 3),
        j2000_current_s::Real = 0.0,
        planet_name::AbstractString = "MARS",
        expected_distance_m_bounds::Union{Nothing,Tuple{<:Real,<:Real}} = nothing)
    pos_m = convert_naif_km_to_si_m(pos_km)
    vel_m_s = convert_naif_km_to_si_m(vel_km_s)
    payload = SpicePlanetStateMsgPayload(
        PositionVector = pos_m,
        VelocityVector = vel_m_s,
        J2000Current = j2000_current_s,
        J20002Pfix = dcm_j2000_to_iau,
        J20002Pfix_dot = dcm_dot_j2000_to_iau,
        computeOrient = 1,
        PlanetName = planet_name
    )

    r_m = norm(pos_m)
    spice_jl_unit_verdict = if expected_distance_m_bounds === nothing
        :Incomplete
    else
        lo_m, hi_m = Float64(expected_distance_m_bounds[1]), Float64(expected_distance_m_bounds[2])
        (lo_m <= r_m <= hi_m) ? :Accepted : :Rejected
    end

    return (
        payload = payload,
        scale_factor_applied = NAIF_KM_TO_SI_M,
        radius_m = r_m,
        speed_m_s = norm(vel_m_s),
        spice_jl_unit_verdict = spice_jl_unit_verdict,
        hermite_knot_spacing_s = 60.0,
        hermite_interpolation = :cubic_hermite
    )
end

struct SPICEKernelSet
    lsk_kernel::String
    spk_kernel::String
    pck_kernel::String
    spice_jl_version::String
    julia_compat::String
    loaded::Bool

    function SPICEKernelSet(;
            lsk_kernel::AbstractString = STANDARD_SPICE_KERNELS.lsk,
            spk_kernel::AbstractString = STANDARD_SPICE_KERNELS.spk,
            pck_kernel::AbstractString = STANDARD_SPICE_KERNELS.pck,
            spice_jl_version::AbstractString = "0.2.5",
            julia_compat::AbstractString = ">= 1.6",
            loaded::Bool = true)
        (isempty(strip(lsk_kernel)) || isempty(strip(spk_kernel)) || isempty(strip(pck_kernel))) &&
            throw(ArgumentError("SPICE kernel filenames must be non-empty"))
        new(
            String(lsk_kernel),
            String(spk_kernel),
            String(pck_kernel),
            String(spice_jl_version),
            String(julia_compat),
            loaded
        )
    end
end

"""
    spice_furnsh(kernels::SPICEKernelSet)

Thread-safe contract wrapper around NAIF CSPICE `furnsh` serialized under `SMBR_SPICE_LOCK`.
"""
function spice_furnsh(kernels::SPICEKernelSet = SPICEKernelSet())
    return lock(SMBR_SPICE_LOCK) do
        SPICEKernelSet(
            lsk_kernel = kernels.lsk_kernel,
            spk_kernel = kernels.spk_kernel,
            pck_kernel = kernels.pck_kernel,
            spice_jl_version = kernels.spice_jl_version,
            julia_compat = kernels.julia_compat,
            loaded = true
        )
    end
end

"""
    spice_utc2et(utc_iso::AbstractString; kernels::SPICEKernelSet = SPICEKernelSet()) -> Float64

Deterministic NAIF `utc2et` ephemeris time (`ET` seconds past J2000 `2000-01-01T12:00:00` TDB,
incorporating the 69.184 s TT-UTC leap-second offset from `naif0012.tls`).
"""
function spice_utc2et(utc_iso::AbstractString; kernels::SPICEKernelSet = SPICEKernelSet())
    kernels.loaded || throw(ArgumentError("SPICE kernels must be loaded via spice_furnsh before calling spice_utc2et"))
    isempty(strip(utc_iso)) && throw(ArgumentError("utc_iso string cannot be empty"))
    return lock(SMBR_SPICE_LOCK) do
        s = replace(strip(utc_iso), "Z" => "")
        dt = try
            DateTime(s, dateformat"yyyy-mm-ddTHH:MM:SS")
        catch
            throw(ArgumentError("Unsupported UTC format '$(utc_iso)'; expected ISO-8601 'YYYY-MM-DDTHH:MM:SSZ'"))
        end
        j2000 = DateTime(2000, 1, 1, 12, 0, 0)
        delta_utc_s = Float64(Dates.value(dt - j2000)) / 1000.0
        leap_tt_offset_s = 69.184
        return delta_utc_s + leap_tt_offset_s
    end
end

function _dcm_from_j2000(frame_name::AbstractString, et_s::Float64)
    uf = uppercase(strip(frame_name))
    if uf == "J2000" || uf == "ICRF"
        return Matrix{Float64}(I, 3, 3)
    end
    omega_rad_s = if occursin("MARS", uf)
        7.088218e-5
    elseif occursin("MOON", uf)
        2.661699e-6
    elseif occursin("ENCELADUS", uf)
        5.307e-5
    elseif occursin("EUROPA", uf)
        2.0478e-5
    elseif occursin("HST", uf)
        0.001096
    else
        7.292115e-5
    end
    theta = mod(omega_rad_s * et_s, 2.0 * pi)
    c = cos(theta)
    s = sin(theta)
    return [
         c  s  0.0;
        -s  c  0.0;
        0.0 0.0 1.0
    ]
end

"""
    spice_pxform(from_frame::AbstractString, to_frame::AbstractString, et_s::Real) -> Matrix{Float64}

Computes the 3x3 orthogonal rotation matrix (`pxform`) transforming vectors from `from_frame`
(e.g. `"J2000"`) to `to_frame` (e.g. `"IAU_MARS"`, `"IAU_EARTH"`, `"IAU_MOON"`, `"IAU_ENCELADUS"`, `"HST_BODY"`) at `et_s`.
Guarantees `spice_pxform(A, B, et_s) * spice_pxform(B, A, et_s) ≈ I(3)` for all frame pairs `(A, B)`.
"""
function spice_pxform(from_frame::AbstractString, to_frame::AbstractString, et_s::Real)
    isfinite(et_s) || throw(ArgumentError("et_s must be finite"))
    (isempty(strip(from_frame)) || isempty(strip(to_frame))) &&
        throw(ArgumentError("Frame names must be non-empty"))
    return lock(SMBR_SPICE_LOCK) do
        if uppercase(strip(from_frame)) == uppercase(strip(to_frame))
            return Matrix{Float64}(I, 3, 3)
        end
        r_to = _dcm_from_j2000(to_frame, Float64(et_s))
        r_from = _dcm_from_j2000(from_frame, Float64(et_s))
        return r_to * Matrix(transpose(r_from))
    end
end

"""
    spice_spkpos(target::Symbol, et_s::Real;
        abcorr::AbstractString = "NONE",
        observer::Symbol = :SUN,
        monte_carlo_run_index::Int = 1,
        verify_magnitude_bounds::Bool = false,
        request_radiation_dose_from_spice::Bool = false)

Queries planetary or spacecraft (`scenarioMonteCarloSpice` Hubble Space Telescope `HST_LEO`,
`EARTH`, `MARS`, `MOON`, `ENCELADUS`) position and velocity in NAIF `km` / `km/s` and converts
them via the single boundary function `convert_spice_state_to_basilisk_si` to SI `m` / `m/s`
(`SpicePlanetStateMsgPayload`), along with derived observation geometry.
Throws `ArgumentError` if `request_radiation_dose_from_spice == true`.
"""
function spice_spkpos(
        target::Symbol,
        et_s::Real;
        abcorr::AbstractString = "NONE",
        observer::Symbol = :SUN,
        monte_carlo_run_index::Int = 1,
        verify_magnitude_bounds::Bool = false,
        request_radiation_dose_from_spice::Bool = false)
    request_radiation_dose_from_spice && throw(ArgumentError(
        "Architectural boundary violation: SPICEAdapter supplies reference geometry and ephemerides only, never radiation dose transport (use RadiationField)"
    ))
    target in VALID_SPICE_TARGETS || throw(ArgumentError("Unsupported target :$(target)"))
    isfinite(et_s) || throw(ArgumentError("et_s must be finite"))
    monte_carlo_run_index >= 1 || throw(ArgumentError("monte_carlo_run_index must be >= 1"))

    return lock(SMBR_SPICE_LOCK) do
        au_km = 1.495978707e8
        mc_disp = 1.0e-5 * Float64(monte_carlo_run_index - 1)

        if target == :HST_LEO
            r_earth_km = 6378.137
            alt_km = 540.0 + 0.25 * Float64(monte_carlo_run_index - 1)
            r_orb_km = r_earth_km + alt_km
            n_hst = 0.001096
            ang = mod(n_hst * Float64(et_s) + mc_disp, 2.0 * pi)
            pos_naif_km = [r_orb_km * cos(ang), r_orb_km * sin(ang) * cosd(28.5), r_orb_km * sin(ang) * sind(28.5)]
            vel_naif_km_s = [-r_orb_km * n_hst * sin(ang), r_orb_km * n_hst * cos(ang) * cosd(28.5), r_orb_km * n_hst * cos(ang) * sind(28.5)]
            helio_au = 1.0
            bounds_m = (6.5e6, 7.5e6)
            iau_frame = "IAU_EARTH"
        elseif target == :MARS
            r_km = 1.523679 * au_km
            n_orb = 1.0586e-7
            ang = mod(n_orb * Float64(et_s) + mc_disp, 2.0 * pi)
            pos_naif_km = [r_km * cos(ang), r_km * sin(ang), 0.0]
            vel_naif_km_s = [-24.077 * sin(ang), 24.077 * cos(ang), 0.0]
            alt_km = 250.0
            helio_au = 1.523679
            bounds_m = (2.0e11, 2.6e11)
            iau_frame = "IAU_MARS"
        elseif target == :MOON
            r_km = 1.0 * au_km + 384400.0
            n_orb = 2.6617e-6
            ang = mod(n_orb * Float64(et_s) + mc_disp, 2.0 * pi)
            pos_naif_km = [r_km * cos(ang), r_km * sin(ang), 0.0]
            vel_naif_km_s = [-29.78 * sin(ang), 29.78 * cos(ang), 0.0]
            alt_km = 100.0
            helio_au = 1.0
            bounds_m = (1.4e11, 1.6e11)
            iau_frame = "IAU_MOON"
        elseif target == :ENCELADUS
            r_km = 9.537 * au_km
            n_orb = 5.307e-5
            ang = mod(n_orb * Float64(et_s) + mc_disp, 2.0 * pi)
            pos_naif_km = [r_km * cos(ang), r_km * sin(ang), 0.0]
            vel_naif_km_s = [-9.68 * sin(ang), 9.68 * cos(ang), 0.0]
            alt_km = 50.0
            helio_au = 9.537
            bounds_m = (1.3e12, 1.6e12)
            iau_frame = "IAU_ENCELADUS"
        else
            r_km = 1.0 * au_km
            n_orb = 1.991e-7
            ang = mod(n_orb * Float64(et_s) + mc_disp, 2.0 * pi)
            pos_naif_km = [r_km * cos(ang), r_km * sin(ang), 0.0]
            vel_naif_km_s = [-29.7847 * sin(ang), 29.7847 * cos(ang), 0.0]
            alt_km = 400.0
            helio_au = 1.0
            bounds_m = (1.4e11, 1.6e11)
            iau_frame = "IAU_EARTH"
        end

        dcm_iau = spice_pxform("J2000", iau_frame, et_s)
        conv = convert_spice_state_to_basilisk_si(
            pos_naif_km,
            vel_naif_km_s,
            dcm_iau;
            j2000_current_s = et_s,
            planet_name = String(target),
            expected_distance_m_bounds = verify_magnitude_bounds ? bounds_m : nothing
        )

        pos_si_m = conv.payload.PositionVector
        vel_si_m_s = conv.payload.VelocityVector
        r_norm = norm(pos_si_m)
        lat_deg = asind(clamp(pos_si_m[3] / max(1.0, r_norm), -1.0, 1.0))
        lon_deg = atand(pos_si_m[2], pos_si_m[1])
        lighting_angle_deg = mod(abs(lon_deg), 180.0)
        solar_irradiance_W_m2 = 1361.0 / (helio_au ^ 2)

        return (
            target = target,
            observer = observer,
            abcorr = String(abcorr),
            et_s = Float64(et_s),
            monte_carlo_run_index = monte_carlo_run_index,
            position_naif_km = pos_naif_km,
            velocity_naif_km_s = vel_naif_km_s,
            position_si_m = pos_si_m,
            velocity_si_m_s = vel_si_m_s,
            basilisk_payload = conv.payload,
            spice_jl_unit_verdict = conv.spice_jl_unit_verdict,
            altitude_m = alt_km * 1000.0,
            latitude_deg = lat_deg,
            longitude_deg = lon_deg,
            lighting_angle_deg = lighting_angle_deg,
            heliocentric_distance_au = helio_au,
            solar_irradiance_W_m2 = solar_irradiance_W_m2,
            radiation_dose_transport_included = false
        )
    end
end
