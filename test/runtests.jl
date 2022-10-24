using FastPKSim
using Test

using ControlSystemsBase

@testset "FastPKSim.jl" begin

    # Test third order model
    # Input data example
    θ_3 = [0.004379828
        0.0036190297
        0.0015290421
        0.0012518701
        8.4043735f-5
        6.7705884]

    # Infusion rates
    u = [0.13888334
        0.13888334
        0.13888334
        0.13888334
        0.13888334
        0.13888334
        0.13888334
        0.13888334
        0.13888334
        0.25
        0.11111667
        0.11111667
        0.11111667
        0.11111667
        0.11111667
        0.11111667
        0.0
        0.0
        0.0
        0.0
        0.0]

    # Bolus doses
    v = [
        130.0
        0.0
        0.0
        0.0
        0.0
        0.0
        0.0
        0.0
        0.0
        0.0
        0.0
        0.0
        0.0
        0.0
        0.0
        0.0
        0.0
        0.0
        0.0
        0.0
        0.0]

    hs = [
        60.0
        60.0
        60.0
        120.0
        300.0
        300.0
        300.0
        600.0
        300.0
        0.08642578
        299.91357
        300.0
        300.0
        600.0
        900.0
        899.8379
        300.1621
        600.0
        900.0
        1800.0
        1800.0]

    youts = [
        2
        3
        4
        5
        6
        7
        8
        9
        12
        13
        14
        15
        16
        18
        19
        20
        21]

    # simulate 
    y = zeros(Float32, length(youts))
    pk3sim!(y, θ_3, u, v, hs, youts)

    # True output
    ytrue = [11.88856
        7.898715
        5.72578
        3.9100208
        3.2476592
        3.2881098
        3.3521297
        3.457583
        3.1101177
        3.0752244
        3.0635502
        3.0553818
        3.064428
        1.3802555
        0.95433724
        0.64292
        0.38235474]

    @test y ≈ ytrue

    # Test if PK3(θ) is correct
    @test all(PK3(θ_3) .≈ (0.1476976506207348, [-6.134363611790634e-5, -0.010055415785499715, -0.0007470542134967702], [-16301.609478739358, -99.4488961303855, -1338.5909374893358], [0.003943514806270174, 0.9436193070471717, 0.05243717814655829]))


    ## Test second order model

    θ_2 = ones(4)
    nt = 100
    time = 0:1:nt
    hs = diff(time)
    u = ones(nt)
    v = ones(nt)
    youts = [1; 10; 13; 14; 15; 16; 20; 21; 22; 78; 79; 82; 84; 87; 90]

    function simulatesecondorder(θ, u, v, h, youts)
        k10, k12, k21, V1 = θ
        A = [-(k10 + k12) k21
            k21 -k21]

        B = [1 / V1; 0]
        C = [1 0]
        D = 0

        PK = ss(A, B, C, D)

        nt = length(u)
        y = zeros(nt,)
        x = [0, 0] # Initial state vector

        for i = 1:nt-1
            if !iszero(v[i])
                x = x + PK.B * v[i] # bolus
                y[i] = x[1]
            end
            if h[i] == 0.0
                y[i+1] = y[i]
                continue
            end
            PK_d = c2d(PK, h[i])
            x = PK_d.A * x + PK_d.B * u[i] # infusion
            y[i+1] = x[1] # infusion affects next sample
        end
        y = y[youts] # observed outputs
        return y
    end

    y = zeros(length(youts))
    pk2sim!(y, θ_2, u, v, hs, youts)

    ytrue = simulatesecondorder(θ_2, u, v, hs, youts)

    @test y ≈ ytrue

    # Test if PK2(θ is correct)
    @test all(PK2(θ_2) .≈ (1.0, [-2.618033988749895, -0.3819660112501051], [-0.38196601125010515, -2.6180339887498953], [0.7236067977499789, 0.27639320225002106]))

end
