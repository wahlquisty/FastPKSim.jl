using FastPKSim
using Test

@testset "FastPKSim.jl" begin
    # Write your tests here.
    @test [2,1] + [1,2] == [3,3]

    @test main() == sqrt.(1:10)

end
