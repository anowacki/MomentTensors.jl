"""
# MomentTensors

The MomentTensors module provides routines for dealing with moment tensors,
such as creating them from double-couple parameters (strike, dip, rake, moment)
or evaluating source radiation patterns from them.

This module uses the Harvard/Global CMT convention for all input/output:
  - x // r         (local radial = upwards)
  - y // θ (theta) (local south = colatitude)
  - z // ϕ (phi)   (local east = (co)longitude)

(Individual routines may use different conventions internally.)

Other conventions are:

* Aki & Richards, Kennett:
  - x // N
  - y // E
  - z // down

MTs are stored internally as immutable types containing only a 6-element array
with the following construction:

    M = [Mrr, Mtt, Mpp, Mrt, Mrp, Mtp]
      = [Mrr, Mθθ, Mϕϕ, Mrθ, Mrϕ, Mθϕ]
    
They can then be accessed using vector index parameters or symbols like so:

    M[:tt] # = 1.5.
    mtp = M[:tp]
    # The following are all true
    M[2,3] == M[:θϕ]
    M[1,2] == Mrθ
    M[1,1] == M.m[1]

All angles are in degrees, always.
"""
module MomentTensors

using LinearAlgebra
using Rotations, StaticArrays

export
    MT,
    amplitude_v_azimuth,
    cmtsolution,
    eps_non_dc,
    m0,
    mw,
    ndk,
    radiation_pattern,
    rotate

"Array to convert from the 6-vector into a 3x3 tensor.  `m.m[_ij2k[i,j]] === m.m[k]`"
const _ij2k = @SMatrix [1 4 5
                        4 2 6
                        5 6 3]

"""
    MT(rr, θθ, ϕϕ, rθ, rϕ, θϕ) -> ::MT
    MT(M::AbstractVector(6)) -> ::MT
    MT(M::AbstractArray(3,3), warn=true) -> ::MT
    MT(strike, dip, rake, M0) -> ::MT

Construct a new `MT` (moment tensor) in the native frame used by `MomentTensors`:
- Radial (`r`) upwards
- Colatitude (`θ` or `t`) southwards
- Longitude (`ϕ` or `p`) eastwards)

Several forms exist to construct a moment tensor:
* Supply individual components as a list of arguments
* Supply a 6-vector
* Specify a strike, dip and rake in degrees, and scalar moment (N.m)
* Give a 3×3 matrix

For the last case, the upper half of the matrix only is used.  By default
symmetry is checked and a warning is issued if the input is not approximately
symmetric.  This can be turned off when `warn` is `false.

One may access the values of a moment tensor `M` in two ways:
1. M[i,j] yields the elements of `M.m` as if they were a two-tensor
2. M[::Symbol] yields the elements by name; see `getindex` for details
"""
struct MT{T<:AbstractFloat}
    m::SVector{6,T}
    MT{T}(m) where T = new{T}(m)
end

MT(args...) = MT{float(promote_type(eltype.(args)...))}(args...)
MT(m::Union{AbstractVector{T}, AbstractArray{T,2}}) where T = MT{float(T)}(m)
MT{T}(rr, tt, pp, rt, rp, tp) where T = MT{T}(SVector{6,T}(rr, tt, pp, rt, rp, tp))
MT{T}(strike, dip, rake, M0) where T = MT{T}(_sdr2mt(strike, dip, rake, M0))
function MT{T}(m::AbstractArray{U,2} where U, warn=true) where T
    size(m) == (3,3) ||
        throw(ArgumentError("2-dimensional array must have dimensions `(3,3)` for a `MT`"))
    warn && !(m[1,2] ≈ m[2,1] && m[1,3] ≈ m[3,1] && m[2,3] ≈ m[3,2]) &&
        @warn("Supplied tensor is not symmetric; using upper elements")
    MT{T}(SVector{6,T}(m[1,1], m[2,2], m[3,3], m[1,2], m[1,3], m[2,3]))
end

"""
    _matrix(m::MT)

Return a 3×3 `SMatrix` for the moment tensor `m`.
"""
_matrix(m::MT) = SMatrix{3,3}(m[1], m[4], m[5],
                              m[4], m[2], m[6],
                              m[5], m[6], m[3])

# Overloaded base operators and constructors
Base.Array(m::MT) = Array(_matrix(m))

"""
    getindex(m::MT, i, j) -> val
    getindex(m::MT, s::Symbol) -> val

Return a component of the moment tensor, `m`.

Two forms are permitted.  Either two indices may be supplied, `i` and `j`, or
a symbol giving the name of the component.  This must be one of:
    :rr, :tt, :pp, :rt, :tr, :rp, :pr, :tp, :pt
"""
Base.getindex(m::MT, s::Symbol) = m.m[_symbol2index(s)]
Base.getindex(m::MT, i) = m.m[i]
Base.getindex(m::MT, i, j) = m.m[_ij2k[i,j]]
Base.getindex(m::MT, i, j, inds...) = m[inds[i],inds[j]][inds...]
Base.size(m::MT) = (3, 3)
Base.isapprox(m1::MT, m2::MT; kwargs...) = isapprox(m1.m, m2.m; kwargs...)

Base.:+(m::MT, a::Real) = MT(m.m .+ a)
Base.:+(a::Real, m::MT) = MT(a .+ m.m)
Base.:+(a::MT, b::MT) = MT(a.m .+ b.m)
Base.:-(m::MT, a::Real) = MT(m.m .- a)
Base.:-(a::Real, m::MT) = MT(a .- m.m)
Base.:-(a::MT, b::MT) = MT(a.m .- b.m)
Base.:-(m::MT) = MT(-m.m)
Base.:*(m::MT, a::Real) = MT(m.m.*a)
Base.:*(a::Real, m::MT) = MT(a.*m.m)
Base.:/(m::MT, a::Real) = MT(m.m./a)
Base.:/(a::Real, m::MT) = MT(a./m.m)

# Routines using or manipulating MTs
"""
    amplitude_v_azimuth(m::MT, inclination, azimuth_range=0:359) -> P, SV, SH, j

Return the P, SV and SH wave amplitudes, `P`, `SV`, `SH` respectively, and source
polarisation, `j`, across the range of azimuth `azimuth_range`.
"""
function amplitude_v_azimuth(m::MT, inclination, azimuth_range=0:359)
    p = Array{Float64}(undef, length(azimuth_range))
    sv = similar(p)
    sh = similar(p)
    j = similar(p)
    for (i, az) in enumerate(azimuth_range)
        p[i], sv[i], sh[i], j[i] = radiation_pattern(m, az, inclination)
    end
    p, sv, sh, j
end

"""
    auxplane(strike, dip, rake) -> strike, dip, rake

Return the strike, dip and rake of the auxiliary plane, given the `strike`,
`dip` and `rake` of the fault plane, all in degrees.
"""
function auxplane(strike, dip, rake)
    # Formula taken from Shearer, Introduction to Seismology
    s1 = deg2rad(strike)
    d1 = deg2rad(dip)
    r1 = deg2rad(rake)

    d2 = acos(sin(r1)*sin(d1))

    sr2 = cos(d1)/sin(d2)
    cr2 = -sin(d1)*cos(r1)/sin(d2)
    r2 = atan(sr2, cr2)

    s12 = cos(r1)/sin(d2)
    c12 = -1/(tan(d1)*tan(d2))
    s2 = s1 - atan(s12, c12)

    strike2 = rad2deg(s2)
    dip2 = rad2deg(d2)
    rake2 = rad2deg(r2)

    if dip2 > 90.0
       strike2 = strike2 + 180.0
       dip2 = 180.0 - dip2
       rake2 = 360.0 - rake2
    end
    rake2 = mod(rake2 + 180.0, 360.0) - 180.0 # In range -180 to 180
    strike2 = mod(strike2, 360.0) # In range 0 to 360
    strike2, dip2, rake2
end

"""
    cmtsolution(str::String) -> m

Return the moment tensor `m` described by the contents of `str`, which are in
the CMTSOLUTION format used by SPECFEM3D(_GLOBE).
"""
function cmtsolution(str::AbstractString)
    lines = _getlines(str)
    length(lines) == 13 || throw(ArgumentError("The supplied string " *
        "does not have 13 lines.  Got:\n" * str))
    # Allow for Fortran double precision scientific notation (1d-3)
    vals = [parse(Float64, replace(split(l)[2], r"[dD]"=>"e")) for l in lines[8:13]]
    MT(vals.*1e-7) # Convert from dyne.cm to N.m
end

"""
    ndk(str::String) -> m, err

Return the moment tensor `m` and its uncertainty `err` given by the contents of
`str`, which are in the ndk format used by the Global CMT project.
"""
function ndk(str::AbstractString)
    lines = _getlines(str)
    length(lines) == 5 || throw(ArgumentError("The supplied string does " *
        "not have 5 lines.  Got:\n" * str))
    tokens = split(lines[4])
    length(tokens) == 13 || throw(ArgumentError("Line 4 of input does not " *
        "have 13 fields.  Got:\n" * lines[4]))
    expo = tryparse(Int, tokens[1])
    expo === nothing && throw(ArgumentError("First token of line 4 " *
        "of input is not an integer.  Got: '$(tokens[1])'"))
    expo = expo - 7 # Convert from dyne.cm to N.m
    m = MT(parse.(Float64, tokens[2:2:end]).*10.0^expo)
    err = MT(parse.(Float64, tokens[3:2:end]).*10.0^expo)
    m, err
end

_getlines(s::String) = split(chomp(s), '\n')

"""
    radiation_pattern(m::MT, azimuth, inclination) -> P, SV, SH, j

Calculate the radiation pattern from the moment tensor `m`, along `azimuth`
measured from local north towards east, and `inclination`, measured away from
downwards towards the radial direction.  Angles in degrees.

### Diagram

    ^ V (SV)
    |
    |     /       View looking along ray
    |_ j /
    | ` /
    | |/
    | /
    |/
    x----------> H (SH)
  
### Returns

* P, SV, SH : Amplitudes of ray in P, SV and SH directionso
* j         : Source polarisation measured from SV towards SH (i.e., like ϕ′ in splitting)
"""
function radiation_pattern(m::MT, azimuth, inclination)
    #=
    Uses formula given in pp. 70 ff. of

        The Seismic Wavefield, Volume 1.  Kennett, B.L.N., Cambridge University Press.

    Kennett uses:

    * x // N
    * y // E
    * z // down

    and:

    - i is angle from z towards x-y plane (incidence angle away from down)
    - phi is angle from x towards y (i.e., azimuth clocwise from N)
    - V is direction upwards (radial) when looking along ray at source
    - H is direction to the right when looking along ray at source:
    =#
    Mxx =  m[:tt]
    Myy =  m[:pp]
    Mzz =  m[:rr]
    Mxy = -m[:tp]
    Mxz =  m[:rt]
    Myz = -m[:rp]

    # Convert to radians
    a = deg2rad(azimuth)
    i = deg2rad(inclination)

    # Some shortcuts
    sini = sin(i)
    cosi = cos(i)
    sinphi = sin(a)
    cosphi = cos(a)
    sintwophi = sin(2a)
    costwophi = cos(2a)
    # Calculate radiation pattern for P, SV and SH
    P = (sini^2)*(Mxx*cosphi^2 + Mxy*sintwophi + Myy*sinphi^2 - Mzz) +
        2*sini*cosi*(Mxz*cosphi + Myz*sinphi) + Mzz
    SV = sini*cosi*(Mxx*cosphi^2 + Mxy*sintwophi + Myy*sinphi^2 - Mzz) +
         cos(2i)*(Mxz*cosphi + Myz*sinphi)
    SH = sini*((Myy - Mxx)*sinphi*cosphi + Mxy*costwophi) +
         cosi*(Myz*cosphi - Mxz*sinphi)
    # Source polarisation in ray frame, measured from upwards towards the right
    j = rad2deg(atan(SV, SH))
    P, SV, SH, j
end

"""
    rotate(m::MT, angle_r, angle_theta, angle_phi) -> ::MT

Return a rotated version of a moment tensor, rotated in turn about the axes
r, theta and phi.  The rotation appears clockwise when looking down each axis
towards the origin.
"""
function rotate(m::MT, r, t, p)
    R = RotZ(deg2rad(p)) * RotY(deg2rad(t)) * RotX(deg2rad(r))
    m′ = R'*Array(m)*R
    MT(m′)
end

"""
    mw(m0) -> mw

Convert a scalar moment (Nm) to a moment magnitude
"""
mw(M0) = 2/3*log10(M0*1e7) - 10.7

"""
    m0(mw) -> m0

Convert a moment magnitude to a scalar moment (Nm)
"""
m0(Mw) = 10^(3/2*(Mw + 10.7))/1e7

"""
    _sdr2mt(strike, dip, rake, m0) -> ::SVector{6}

Return a moment tensor defined by the `strike`, `dip` and `rake` of the fault
plane, and a scalar moment `m0`.  Use `m0(mw)` if you want to specify the
moment mangitude instead.

Convention:
    * Strike is measured such that the downdip direction is 90° greater in azimuth
      than the strike.
    * Dip is measured between 0 and 90°.
    * Rake is measured from the strike direction, positive anticlockwise on the
      footwall, such that positive rakes imply a thrust component, and negative
      rakes a normal component.  Rakes should be in the range -180° to 180°

All angles are in degrees.
"""
function _sdr2mt(strike, dip, rake, M0)
    0.0 .<= dip .<= 90.0 || throw(ArgumentError("`dip` must be between 0° and 90°"))
    s, d, r = deg2rad(strike), deg2rad(dip), deg2rad(rake)
    tt = -M0*(sin(d)*cos(r)*sin(2s) + sin(2d)*sin(r)*sin(s)^2)
    tp = -M0*(sin(d)*cos(r)*cos(2s) + sin(2d)*sin(r)*sin(2s)/2)
    rt = -M0*(cos(d)*cos(r)*cos(s)  + cos(2d)*sin(r)*sin(s))
    pp =  M0*(sin(d)*cos(r)*sin(2s) - sin(2d)*sin(r)*cos(s)^2)
    rp =  M0*(cos(d)*cos(r)*sin(s)  - cos(2d)*sin(r)*cos(s))
    rr =  M0*sin(2d)*sin(r)
    @SVector [rr, tt, pp, rt, rp, tp]
end

"""
    _symbol2index(s::Symbol) -> k

Return the index in the range `1 <= k <= 6` with which to access the `m` field
of an `MT` type.

The following are permitted; case and order do not matter:
    :rr, :tt, :θθ, :pp, :ϕϕ, :rt, :rθ, :rp, :rϕ, :tp, :θϕ
"""
function _symbol2index(s::Symbol)
    if     s == :rr                  1
    elseif s in (:tt, :θθ)           2
    elseif s in (:pp, :ϕϕ)           3
    elseif s in (:rt, :tr, :rθ, :θr) 4
    elseif s in (:rp, :pr, :rϕ, :ϕr) 5
    elseif s in (:tp, :pt, :θϕ, :ϕθ) 6
    else throw(ArgumentError("'$s' is not a valid MT component name"))
    end
end

"""
    eps_non_dc(m::MT) -> ϵ

Return the deviatoric non-double-couple component `ϵ` of the moment tensor `m`,
defined by Giardini (1983) as

ϵ = –λ₂/max(|λ₁|, |λ₃|),

where the eigenvalues of the moment tensor are λ₁ ≥ λ₂ ≥ λ₃.

## References
- Giardini, D., 1983. Regional deviation of earthquake source mechanisms from the
  'double couple' model. In: H. Kanamori and E. Boschi (Editors), Earthquakes:
  Observation North-Holland, Amsterdam, Theory and Interpretation, pp. 345-353.
"""
function eps_non_dc(m::MT)
    λ = sort(eigen(Symmetric(_matrix(m))).values)
    -λ[2]/max(abs(λ[1]), abs(λ[3]))
end

# Printing to the REPL looks like a 3×3 matrix
function Base.show(io::IO, mime::MIME"text/plain", m::MT{T}) where T
    println(io, "MomentTensors.MT{$T}:")
    io = IOContext(io, :compact=>true)
    Base.print_array(io, Array(m))
end

end # module
