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

using Rotations

import Base: +, -, *, /, Array, endof, getindex, ndims, size

export
    MT,
    m0,
    mw,
    radiation_pattern,
    rotate

"Array to convert from the 6-vector into a 3x3 tensor.  `m.m[_ij2k[i,j]] === m.m[k]`"
const _ij2k = [1 4 5
               4 2 6
               5 6 3]

"""
    MT(rr, θθ, ϕϕ, rθ, rϕ, θϕ) -> ::MT
    MT(M::Vector(6)) -> ::MT
    MT(M::Array(3,3)) -> ::MT
    MT(strike, dip, rake, M0) -> ::MT

Construct a new `MT` (moment tensor) in the native frame used by `MomentTensors`:
- Radial (`r`) upwards
- Colatitude (`θ` or `t`) southwards
- Longitude (`ϕ` or `p`) eastwards)

Several forms exist to construct a moment tensor:
* Supply individual components as a list of arguments
* Supply a 6-vector
* Give a 3×3 matrix
* Specify a strike, dip and rake in degrees, and scalar moment (N.m)

The `MT` type holds one field, `m`, as a length-6 vector.

One may access the values of a moment tensor `M` in two ways (beyond directly
accessing the field `M.m`):
1. M[i,j] yields the elements of `M.m` as if they were a two-tensor
2. M[::Symbol] yields the elements by name; see `getindex` for details
"""
immutable MT
    m :: Vector{Float64}
    MT{T}(m::Vector{T}) = new(m)
    MT(rr, tt, pp, rt, rp, tp) = new([rr, tt, pp, rt, rp, tp])
    MT{T}(m::Array{T,2}) = size(m) == (3,3) &&
        new([m[1,1], m[2,2], m[3,3], m[1,2], m[1,3], m[2,3]]) ||
        error("2-dimensional array must have dimensions `(3,3)` for a `MT`")
    MT(strike, dip, rake, M0) = new(_sdr2mt(strike, dip, rake, M0))
end

# Overloaded base operators and constructors
Array(m::MT) = [m[i,j] for i in 1:3, j in 1:3]

"""
    getindex(m::MT, i, j) -> val
    getindex(m::MT, s::Symbol) -> val

Return a component of the moment tensor, `m`.

Two forms are permitted.  Either two indices may be supplied, `i` and `j`, or
a symbol giving the name of the component.  This must be one of:
    :rr, :tt, :pp, :rt, :tr, :rp, :pr, :tp, :pt
"""
getindex(m::MT, s::Symbol) = m.m[_symbol2index(s)]
getindex(m::MT, i::Integer, j::Integer) = m.m[_ij2k[i,j]]
getindex(m::MT, inds...) = getindex(Array(m), inds...)
ndims(m::MT) = 2
size(m::MT) = (3, 3)
function size(m::MT, dim)
    if dim < 1
        error("arraysize: dimension out of range")
    elseif dim < 3
        3
    else
        1
    end
end
endof(m::MT) = endof(Array(m))

+(m::MT, a::Real) = MT(m.m .+ a)
+(a::Real, m::MT) = m + a
-(m::MT, a::Real) = MT(m.m .- a)
-(a::Real, m::MT) = MT(a .- m.m)
*(m::MT, a::Real) = MT(m.m.*a)
*(a::Real, m::MT) = m*a
/(m::MT, a::Real) = MT(m.m./a)
/(a::Real, m::MT) = MT(a./m.m)

# Routines using or manipulating MTs
"""
    auxplane(strike, dip, rake) -> strike, dip, rake

Return the strike, dip and rake of the auxiliary plane, given the `strike`,
`dip` and `rake` of the fault plane, all in degrees.
"""
function auxplane(strike, dip, rake)
    # Formula taken from Shearer, Introduction to Seismology
    s1 = deg2rad(strike1)
    d1 = deg2rad(dip1)
    r1 = deg2rad(rake1)

    d2 = acos(sin(r1)*sin(d1))

    sr2 = cos(d1)/sin(d2)
    cr2 = -sin(d1)*cos(r1)/sin(d2)
    r2 = atan2(sr2,cr2)

    s12 = cos(r1)/sin(d2)
    c12 = -1/(tan(d1)*tan(d2))
    s2 = s1 - atan2(s12, c12)

    strike2 = rad2deg(s2)
    dip2 = rad2deg(d2)
    rake2 = rad2deg(r2)

    if dip2 > 90.
       strike2 = strike2 + 180.
       dip2 = 180. - dip2
       rake2 = 360. - rake2
    end
    rake2 = mod(rake2 + 180., 360.) - 180. # In range -180 to 180
    strike2 = mod(strike2, 360.) # In range 0 to 360
    strike2, dip2, rake2
end

"""
    radiation_pattern(m::MT, azimuth, inclination) -> P, SV, SH, j

Calculate the radiation pattern from the moment tensor `m`, along `azimuth`
measured from local north towards east, and `inclination`, measured away from
downwards towards the radial direction.  Angles in degrees.

### Diagram

    ^ V (SV)
    |
    | j   /       View looking along ray
    |-   /
    | \ /
    |  /
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
 	sintwophi = sin(2*a)
 	costwophi = cos(2*a)
 	# Calculate radiation pattern for P, SV and SH
 	P = (sini^2)*(Mxx*cosphi^2 + Mxy*sintwophi + Myy*sinphi^2 - Mzz) +
 		2*sini*cosi*(Mxz*cosphi + Myz*sinphi) + Mzz
 	SV = sini*cosi*(Mxx*cosphi^2 + Mxy*sintwophi + Myy*sinphi^2 - Mzz) +
 		 cos(2*i)*(Mxz*cosphi + Myz*sinphi)
 	SH = sini*((Myy - Mxx)*sinphi*cosphi + Mxy*costwophi) +
 		 cosi*(Myz*cosphi - Mxz*sinphi)
 	# Source polarisation in ray frame, measured from upwards towards the right
 	j = rad2deg(atan2(SV,SH))
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
    MT(Array(m′))
end

"    mw(m0) -> mw\n\nConvert a scalar moment (Nm) to a moment magnitude"
mw(M0) = (2/3)*(log10(M0) - 16.1)
"    m0(mw) -> m0\n\nConvert a moment magnitude to a scalar moment (Nm)"
m0(Mw) = 10^((3/2)*Mw + 16.1)

"""
    _sdr2mt(strike, dip, rake, m0) -> ::MT

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
    0. .<= dip .<= 90. || error("`dip` must be between 0° and 90°")
    s, d, r = deg2rad(strike), deg2rad(dip), deg2rad(rake)
    tt = -M0*(sin(d)*cos(r)*sin(2s) + sin(2d)*sin(r)*sin(s)^2)
    tp = -M0*(sin(d)*cos(r)*cos(2s) + sin(2d)*sin(r)*sin(2s)/2)
    rt = -M0*(cos(d)*cos(r)*cos(s)  + cos(2d)*sin(r)*sin(s))
    pp =  M0*(sin(d)*cos(r)*sin(2s) - sin(2d)*sin(r)*cos(s)^2)
    rp =  M0*(cos(d)*cos(r)*sin(s)  - cos(2d)*sin(r)*cos(s))
    rr =  M0*sin(2d)*sin(r)
    [rr, tt, pp, rt, rp, tp]
end

"""
    _symbol2index(s::Symbol) -> k

Return the index in the range `1 <= k <= 6` with which to access the `m` field
of an `MT` type.

The following are permitted; case and order do not matter:
    :rr, :tt, :θθ, :pp, :ϕϕ, :rt, :rθ, :rp, :rϕ, :tp, :θϕ
"""
function _symbol2index(s::Symbol)
    s = Symbol(lowercase(string(s)))
    if     s == :rr                  1
    elseif s in (:tt, :θθ)           2
    elseif s in (:pp, :ϕϕ)           3
    elseif s in (:rt, :tr, :rθ, :θr) 4
    elseif s in (:rp, :pr, :rϕ, :ϕr) 5
    elseif s in (:tp, :pt, :θϕ, :ϕθ) 6
    else error("'$s' is not a valid MT component name")
    end
end

end # module
