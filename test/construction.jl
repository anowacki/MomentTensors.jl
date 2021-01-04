using MomentTensors, Test, StaticArrays

@testset "Construction" begin
    @testset "Direct" begin
        @test MT(1, 2, 3, 4, 5, 6) isa MT{Float64}
        @test MT{Float32}(1, 2, 3, 4, 5, 6) isa MT{Float32}
        @testset "$T" for T in (Float32, Float64)
            @test eltype(MT{T}(1, 2, 3, 4, 5, 6)) == T
            @test eltype(MT{Complex{T}}(1, 2, 3, 4, 5, 6)) == Complex{T}
        end
        @test MT(1, 2, 3, 4, 5, 6) == MT(@SVector[1, 2, 3, 4, 5, 6])
        @test MT(Int32[1, 2, 3, 4, 5, 6]) ==
            MT{float(Int32)}(@SVector[1.f0, 2.f0, 3.f0, 4.f0, 5.f0, 6.f0])
        @test MT(Float32[1, 2, 3, 4, 5, 6]) == MT{Float32}(1, 2, 3, 4, 5, 6)
        @test MT([11 12 13; 12 22 23; 13 23 33]) == MT(11, 22, 33, 12, 13, 23)
        @test MT(11, 22, 33, 12, 13, 23) == MT([11, 22, 33, 12, 13, 23])
        @test_throws ArgumentError MT([1 2 3; 4 5 6])
        # Warn for asymmetric matrices
        @test_logs (:warn,
            "Supplied tensor is not symmetric; using upper elements") MT([1 2 3
                                                                          4 5 6
                                                                          7 8 9])
        # Do not warn with `warn=false`
        @test_logs MT([1 2 3; 4 5 6; 7 8 9], false)
    end

    @testset "Strike-dip-rake" begin
        @test MT(30, 60, -45, 7e15) ≈
            MT(1e16.*(-0.429, -0.264, 0.693, -0.338, -0.091, -0.029)...) atol=0.01e16
    end
end