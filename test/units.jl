# MTs with eltypes which have units from the Unitful package

using Test
using MomentTensors
using Unitful

@testset "Units" begin
    m = Unitful.m
    N = Unitful.N
    Nm = N*m

    @testset "Construction" begin
        @test MT((1, 2, 3, 4, 5, 6).*u"N*m") == MT(1Nm, 2Nm, 3Nm, 4Nm, 5Nm, 6Nm) == MT((1:6)*Nm)
        @test MT(1:6) != MT((1:6)*Nm)
        let mt = MT((1:6)*Nm)
            @test eltype(mt) == typeof(float(1Nm))
        end
    end

    @testset "getindex/setindex!" begin
        let mt = MT((1:6)*Nm)
            @test mt[1,1] == 1Nm
        end
    end

    @testset "Arithmetic" begin
        @test MT(ones(6).*u"N*m/s") + MT(2*ones(6).*u"N*m/s") == MT(3*ones(6).*u"N*m/s")
        let m1 = MT((11, 22, 33, 12, 13, 23).*Nm),
            m2 = MT((-11, -22, -33, -12, -13, -23).*Nm)
            @test m1 + m2 == MT((0, 0, 0, 0, 0, 0).*Nm)
            @test m1 == -m2
            @test -m1 == m2
            @test m1 + 1Nm == MT((12, 23, 34, 13, 14, 24).*Nm)
            @test m1 - 1Nm == MT((10, 21, 32, 11, 12, 22).*Nm)
            @test m1 + 3Nm == 3Nm + m1
            @test m1 - 4Nm == -(4Nm - m1)
            @test 2m1 == MT((22, 44, 66, 24, 26, 46).*Nm)
            @test 2m1 == m1*2
            @test 1/m1 == MT(1 ./ ((11, 22, 33, 12, 13, 23).*Nm)...)
            @test m1/2 == MT((5.5, 11, 16.5, 6, 6.5, 11.5).*Nm)
        end
    end

    @testset "Conversion" begin
        @test Array(MT((1:6)*Nm)) == Float64[1 4 5; 4 2 6; 5 6 3].*Nm
    end

    @testset "Decomposition" begin
        let m = MT([ 3 -5 10
                    -5  1  4
                    10  4  1].*Nm)
            eltp = typeof(float(1Nm))
            d = decompose(m)

            @testset "Unitful MT $f" for f in (:iso, :dev, :dc, :clvd)
                @test eltype(getfield(d, f)) == eltp
            end
            @testset "Unitful scalar $f" for f in (:iso_m0, :dev_m0, :m0)
                @test typeof(getfield(d, f)) == eltp
            end
            @testset "Dimensionless scalar $f" for f in (:prop_iso, :prop_dev, :prop_dc, :prop_clvd)
                @test typeof(getfield(d, f)) == typeof(ustrip(m[1,1]))
            end
        end
    end

    @testset "amplitude_v_azimuth" begin
        eltp = typeof(Float32(1N))
        p, sv, sh, j = amplitude_v_azimuth(MT(rand(Float32, 6).*N), 0, 0)
        for v in (p, sv, sh)
            @test eltype(v) == eltp
        end
        @test eltype(j) == Float32
    end

    @testset "radiation_pattern" begin
        eltp = typeof(Float32(1Nm))
        p, sv, sh, j = radiation_pattern(MT(rand(Float32, 6).*Nm), 0, 0)
        for v in (p, sv, sh)
            @test v isa eltp
        end
        @test j isa Float32
    end

    @testset "rotate" begin
        @test rotate(MT(rand(Float32, 6).*Nm), 1, 2, 3) isa MT{typeof(1f0*Nm)}
    end

    @testset "eps_non_dc" begin
        @test eps_non_dc(MT(rand(Float32, 6).*Nm)) isa typeof(ustrip(1f0*Nm))
    end
end
