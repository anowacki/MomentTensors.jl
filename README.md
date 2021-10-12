# MomentTensors.jl

[![Build Status](https://github.com/anowacki/MomentTensors.jl/workflows/CI/badge.svg)](https://github.com/anowacki/MomentTensors.jl/actions)
[![Coverage status](https://codecov.io/gh/anowacki/MomentTensors.jl/branch/master/graph/badge.svg?token=knbujQ671A)](https://codecov.io/gh/anowacki/MomentTensors.jl)

## What is MomentTensors.jl?
A [Julia](http://julialang.org) package for dealing with [seismic moment
tensors](https://earthquake.usgs.gov/learn/glossary/?term=moment%20tensor).

It is currently very limited, and useful for three main things:

1. Calculating the radiation pattern for a moment tensor;
2. Decomposing tensors into their istropic, double-couple and CLVD
   components; and
3. Rotating moment tensors.

I wrote this module because these are by far the most common things I need to do
in my day-to-day global seismologist life.

## What it isn't
Currently, no plotting is performed, nor conversion between conventions.

The module internally assumes the Harvard/Global CMT convention (see the module
interactive help for details), though nothing is stopping you using another
convention with it as long as you remember which indices correspond to which
directions.

Input/output is limited to reading the CMTSOLUTION format used by
[SPECFEM3D](https://github.com/geodynamics/specfem3d) and the
[NDK](http://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/allorder.ndk_explained)
format used by the [GlobalCMT project](https://globalcmt.org).

## How to install
MomentTensors.jl can be added to your Julia install like so:

```julia
julia> import Pkg; Pkg.pkg"add https://github.com/anowacki/MomentTensors.jl"
```


## How to use
### MT type
MomentTensors.jl represents moment tensors using the `MT` type.

#### Construction
They can be constructed by:
- passing a set of six moment tensor components, in the order
  _M<sub>rr</sub>, M<sub>θθ</sub>, M<sub>φφ</sub>,
  M<sub>rθ</sub>, M<sub>rφ</sub>, M<sub>θφ</sub>_),
- a vector of length six containing these components,
- a 3 × 3 matrix _M<sub>ij</sub>_ with
  _i,j_ ∈ {_r,θ,φ_}, or
- the strike, dip and rake (in Aki & Richards convention) plus a moment
  in N.m.

The element type of an `MT{T} where T` instance is determined
automatically from the values supplied, or one can specify a
desired element type explicitly with e.g. `MT{Float32}(1, 2, 3, 4, 5, 6)`.

See the docstring for `MT` for more details of construction.

#### Indexing
To retrieve an individual component of the `MT` `m`, you can:
- access the _i,j_ component with `m[i,j]` where _i,j_ ∈ {1,2,3};
- get the named components with `m[:rθ]` or equivalently `m[:rt]`
  (see the docstring for `getindex`); and
- get the _I_<sup>th</sup> component of the six-element vector with
  `m[I]`.


## Exported functions

- `MT`: Construct a new moment tensor.
- `amplitude_v_azimuth`: Compute the P, SV and SH amplitudes, and polarisation angle,
  for a particular takeoff angle at a range of azimuths.
- `cmtsolution`: Construct a new moment tensor from a string in the SPECFEM3D 'CMTSOLUTION'
   format
- `decompose`: Decompose a moment tensor into its isotropic, double-couple
  and CLVD components, and report the relative proportion of the isotropic,
  deviatoric and double-couple parts, plus their associated moments.
- `eps_non_dc`: Calculate the non-double-couple component of an MT as
  defined by Giardini.
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
      

  One may access the values of a moment tensor M in two ways:

    1.   M[i,j] yields the elements of M.m as if they were a two-tensor
      
    2.   M[::Symbol] yields the elements by name; see getindex(::MT) for details
```

## Contributing
If you find a bug with MomentTensors or have suggestions for improvement,
please
[open an issue](https://github.com/anowacki/MomentTensors.jl/issues/new/choose)
giving as much information as possible on how to reproduce the bug or problem.

[Pull requests](https://github.com/anowacki/MomentTensors.jl/compare)
to add new features are welcome and will be seriously
considered.  Please note that the package aims to be lightweight and
rely on few external dependencies.
