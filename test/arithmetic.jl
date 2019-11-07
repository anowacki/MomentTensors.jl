using MomentTensors, Test

@testset "Operations" begin
    let m1 = MT(11, 22, 33, 12, 13, 23), m2 = MT(-11, -22, -33, -12, -13, -23)
        @test m1 + m2 == MT(0, 0, 0, 0, 0, 0)
        @test m1 == -m2
        @test -m1 == m2
        @test m1 + 1 == MT(12, 23, 34, 13, 14, 24)
        @test m1 - 1 == MT(10, 21, 32, 11, 12, 22)
        @test m1 + 3 == 3 + m1
        @test m1 - 4 == -(4 - m1)
        @test 2m1 == MT(22, 44, 66, 24, 26, 46)
        @test 2m1 == m1*2
        @test 1/m1 == MT(1 ./ (11, 22, 33, 12, 13, 23)...)
        @test m1/2 == MT(5.5, 11, 16.5, 6, 6.5, 11.5)
    end
end
