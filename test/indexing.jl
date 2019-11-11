using MomentTensors, Test

@testset "Indexing" begin
    let v = [11, 22, 33, 12, 13, 23], m = MT(v),
            I = [1 4 5; 4 2 6; 5 6 3], a = Array(m)
        for i in 1:6
            @test m[i] == v[i]
        end
        for i in 1:3, j in i:3
            # Symmetry
            @test m[i,j] == m[j,i]
            @test m[i,j] == a[i,j]
            @test m[i,j] == a[j,i]
            # Correct value
            @test m[i,j] == 10i + j
            @test a[i,j] == 10i + j
            # Correct vector index
            @test m[i,j] == m[I[i,j]]
            @test m[i,j] == m.m[I[i,j]]
        end
        for (i, s1) in enumerate((:r, :t, :p)), (j, s2) in enumerate((:r, :t, :p))
            @test m[i,j] == m[Symbol(s1, s2)]
        end
        for (i, s1) in enumerate((:r, :θ, :ϕ)), (j, s2) in enumerate((:r, :θ, :ϕ))
            @test m[i,j] == m[Symbol(s1, s2)]
        end
        @test size(m) == (3, 3)
        @test_throws ArgumentError m[:xy]
        @test_throws BoundsError m[7]
        @test_throws BoundsError m[1,4]
    end
end
