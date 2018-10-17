# MomentTensors.jl

## What is MomentTensors.jl?
A [Julia](http://julialang.org) package for dealing with [seismic moment
tensors](https://earthquake.usgs.gov/learn/glossary/?term=moment%20tensor).

It is currently very limited, and useful for two main things:

1. Calculating the radiation pattern for a moment tensor; and
2. Rotating moment tensors.

I wrote this module because these are by far the most common things I need to do
in my day-to-day global seismologist life.

## What it isn't
Currently, no plotting is performed, nor conversion between conventions.

The module internally assumes the Harvard/Global CMT convention (see the module
interactive help for details), though nothing is stopping you using another
convention with it as long as you remember which indices correspond to which
directions.

No input/output is performed beyond basic Julia IO.

## How to install
Although not registered as an official package, MomentTensors.jl can be added
to your Julia install like so:

```julia
julia> import Pkg; Pkg.add("https://github.com/anowacki/MomentTensors.jl")
```

(On Julia versions less than v0.7, do:
`Pkg.clone("https://github.com/anowacki/MomentTensors.jl"); Pkg.checkout("MomentTensors", "julia0.6")`.)

You then need only do

```julia
using MomentTensors
```

and if that works, you're ready to go.


## How to use
### MT type
MomentTensors.jl represents moment tensors using the `MT` type.  This is an
immutable type with one field (`m`) which is a length-6 `StaticArrays.SVector`.

Moment tensors can be specified by one of several ways; see the example help
output in the 'Getting help' section below for the possible invocations.

## Exported functions

- `MT`: Construct a new moment tensor.
- `amplitude_v_azimuth`: Compute the P, SV and SH amplitudes, and polarisation angle,
  for a particular takeoff angle at a range of azimuths.
- `cmtsolution`: Construct a new moment tensor from a string in the SPECFEM3D 'CMTSOLUTION'
   format
- `eps_non_dc`: Calculate the non-double-couple component of an MT.
- `m0`: Return the scalar moment, given a moment magnitude.
- `mw`: Return the moment magnitude, given a scalar moment.
- `ndk`: Construct a new moment tensor from a string in the 'NDK' format used by
  the Global CMT project.
- `radiation_pattern`: Compute the P, SV and SH amplitude, and S polarisation angle,
  along a specific takeoff angle and azimuth.
- `rotate`: Rotate an MT.


## Getting help
Functions are documented, so at the REPL type `?` to get a `help?>` prompt,
and type the name of the function:

```julia
help?> MT
search: MT mtime SymTridiagonal Meta Method match Matrix mktemp methods matchall

  MT(rr, θθ, ϕϕ, rθ, rϕ, θϕ) -> ::MT
  MT(M::Vector(6)) -> ::MT
  MT(M::Array(3,3)) -> ::MT
  MT(strike, dip, rake, M0) -> ::MT

  Construct a new MT (moment tensor) in the native frame used by MomentTensors:

    •    Radial (r) upwards
      
    •    Colatitude (θ or t) southwards
      
    •    Longitude (ϕ or p) eastwards)
      

  Several forms exist to construct a moment tensor:

    •    Supply individual components as a list of arguments
      
    •    Supply a 6-vector
      
    •    Give a 3×3 matrix
      
    •    Specify a strike, dip and rake in degrees, and scalar moment (N.m)
      

  The MT type holds one field, m, as a length-6 vector.

  One may access the values of a moment tensor M in two ways (beyond directly accessing
  the field M.m):

    1.   M[i,j] yields the elements of M.m as if they were a two-tensor
      
    2.   M[::Symbol] yields the elements by name; see getindex for details
```
