using MomentTensors, Test

@testset "I/O" begin
    @testset "GlobalCMT" begin
        let (m, e) = ndk("""
            PDE  2000/01/01 05:24:35.3  36.87   69.95  54.3 5.2 4.5 HINDU KUSH REGION       
            B010100A         B: 13   21  45 S:  0    0   0 M:  0    0   0 CMT: 1 BOXHD:  1.0
            CENTROID:      4.0 1.0  37.15 0.13   69.72 0.14  84.6  7.2 FREE O-00000000000000
            23  2.495 0.595 -1.952 0.880 -0.544 0.633 -0.154 0.774  0.416 0.840 -5.825 1.128
            V10   4.700 11 228   2.420 79  43  -7.120  1 138   5.910 273 82  173   4 83    8
            """)
            @test m ≈ MT(1e16.*[2.495, -1.952, -0.544, -0.154, 0.416, -5.825]) atol=1e13
            @test e ≈ MT(1e16.*[0.595,  0.880,  0.633,  0.774, 0.840,  1.128]) atol=1e13
        end
        # Wrong number of columns on fourth line
        @test_throws ArgumentError ndk("""
            PDE  2000/01/01 05:24:35.3  36.87   69.95  54.3 5.2 4.5 HINDU KUSH REGION       
            B010100A         B: 13   21  45 S:  0    0   0 M:  0    0   0 CMT: 1 BOXHD:  1.0
            CENTROID:      4.0 1.0  37.15 0.13   69.72 0.14  84.6  7.2 FREE O-00000000000000
            23  2.495 0.595 -1.952 0.880 -0.544 0.633 -0.154 0.774  0.416 0.840 -5.825
            V10   4.700 11 228   2.420 79  43  -7.120  1 138   5.910 273 82  173   4 83    8
            """)
        # Incorrect number of lines
        @test_throws ArgumentError ndk("""
            B010100A         B: 13   21  45 S:  0    0   0 M:  0    0   0 CMT: 1 BOXHD:  1.0
            CENTROID:      4.0 1.0  37.15 0.13   69.72 0.14  84.6  7.2 FREE O-00000000000000
            23  2.495 0.595 -1.952 0.880 -0.544 0.633 -0.154 0.774  0.416 0.840 -5.825 1.128
            V10   4.700 11 228   2.420 79  43  -7.120  1 138   5.910 273 82  173   4 83    8
            """)
        # Exponent is not an integer
        @test_throws ArgumentError ndk("""
            PDE  2000/01/01 05:24:35.3  36.87   69.95  54.3 5.2 4.5 HINDU KUSH REGION       
            B010100A         B: 13   21  45 S:  0    0   0 M:  0    0   0 CMT: 1 BOXHD:  1.0
            CENTROID:      4.0 1.0  37.15 0.13   69.72 0.14  84.6  7.2 FREE O-00000000000000
            23.0  2.495 0.595 -1.952 0.880 -0.544 0.633 -0.154 0.774  0.416 0.840 -5.825 1.128
            V10   4.700 11 228   2.420 79  43  -7.120  1 138   5.910 273 82  173   4 83    8
            """)
    end

    @testset "CMTSOLUTION" begin
        @test cmtsolution("""
            PDE 1994  6  9  0 33 16.40 -13.8300  -67.5600 637.0 6.9 6.8 NORTHERN BOLIVIA              
            event name:     060994A        
            time shift:     29.0000
            half duration:  20.0000
            latitude:      -13.8200
            longitude:     -67.2500
            depth:         647.1000
            Mrr:      -7.590000d+27
            Mtt:       7.750000d+27
            Mpp:      -1.600000e+26
            Mrt:      -2.503000e+28
            Mrp:       4.200000e+26
            Mtp:      -2.480000e+27
            """) ≈ MT(-7.59e20, 7.75e20, -1.6e19, -2.503e21, 4.2e19, -2.48e20) atol=1e14
        # Missing a line
        @test_throws ArgumentError cmtsolution("""
            PDE 1994  6  9  0 33 16.40 -13.8300  -67.5600 637.0 6.9 6.8 NORTHERN BOLIVIA              
            event name:     060994A        
            time shift:     29.0000
            half duration:  20.0000
            latitude:      -13.8200
            longitude:     -67.2500
            depth:         647.1000
            Mrr:      -7.590000e+27
            Mtt:       7.750000e+27
            Mpp:      -1.600000e+26
            Mrt:      -2.503000e+28
            Mrp:       4.200000e+26
            """)
    end
end