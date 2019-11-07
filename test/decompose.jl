using MomentTensors, Test

# Comparison with example at http://www.larskrieger.de/mopad/#getall
@testset "Decomposition" begin
    let m = MT([ 3 -5 10
                -5  1  4
                10  4  2])
        d = decompose(m)
        @test d.iso == MT(2, 2, 2, 0, 0, 0)
        @test d.prop_iso ≈ 0.13 atol=0.01
        @test d.dev == MT([ 1 -5 10
                           -5 -1  4
                           10  4  0])
        @test d.prop_dev ≈ 0.87 atol=0.01
        @test d.dc ≈ MT([1.36 -2.97  7.30
                        -2.97 -1.77  1.95
                         7.30  1.95  0.41]) atol=0.01
        @test d.prop_dc ≈ 0.55 atol=0.01
        @test d.clvd ≈ MT([-0.36 -2.03  2.70
                           -2.03  0.77  2.05
                            2.70  2.05 -0.41]) atol=0.01
        @test d.prop_clvd ≈ 0.32 atol=0.01
        @test d.m0 ≈ 14.9031089939
        @test mw(d.m0) ≈ -5.25114874821
    end
end
