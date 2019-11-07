using MomentTensors, Test

@testset "Planes" begin
    @test all(isapprox.(MomentTensors.auxplane(30, 60, -45),
            (147.0, 52.0, -141.0), atol=1))
end