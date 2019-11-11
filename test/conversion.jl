using MomentTensors, Test

@testset "Conversion" begin
    @test m0(mw(1e18)) ≈ 1e18
    @test mw(m0(5)) ≈ 5
end
