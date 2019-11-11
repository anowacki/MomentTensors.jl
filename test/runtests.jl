using MomentTensors, Test

@testset "All tests" begin
    include("construction.jl")
    include("indexing.jl")
    include("arithmetic.jl")
    include("conversion.jl")
    include("planes.jl")
    include("io.jl")
    include("decompose.jl")
end
