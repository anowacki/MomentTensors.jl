using MomentTensors, Test, StaticArrays

# Comparison with example at http://www.larskrieger.de/mopad/#getall
@testset "Decomposition" begin
    let m = MT([ 3 -5 10
                -5  1  4
                10  4  2])
        d = decompose(m)
        @testset "Has field $f" for f in (
                :iso, :dev, :dc, :clvd, :iso_m0, :dev_m0, :prop_iso, :prop_dev,
                :prop_dc, :prop_clvd, :m0)
            @test f in fieldnames(typeof(d))
        end
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
        @test eps_non_dc(m) ≈ -0.34250995536253115
    end

    @testset "_sortperm_abs_3" begin
        @test MomentTensors._sortperm_abs_3(SVector(1.0, 2.0, 3.0)) == SVector(1, 2, 3)
        @test MomentTensors._sortperm_abs_3(SVector(-3.0, -2.0, 1)) == SVector(3, 2, 1)
    end
end
