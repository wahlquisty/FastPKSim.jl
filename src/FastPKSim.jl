module FastPKSim

export main, pk3sim!, pk2sim!, PK2, PK3

using LinearAlgebra,StaticArrays, SLEEFPirates

include("simulation.jl") # functions to simulate state and compute output
include("pkmodels.jl") # calculations of λ, R from parameter vector θ

end