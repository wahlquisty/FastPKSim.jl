module FastPKSim

export main, pk3sim!

using LinearAlgebra,StaticArrays, SLEEFPirates

include("simulation.jl") # functions to simulate state and compute output
include("pkmodels.jl") # calculations of λ, R from parameter vector θ

end