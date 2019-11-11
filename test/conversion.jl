using MomentTensors, Test

@testset "Conversion" begin
    @test m0(mw(1e18)) ≈ 1e18
    @test mw(m0(5)) ≈ 5
    @test Array(MT{Float16}(0, 45, 90, 1)) isa Array{Float16,2}
end
