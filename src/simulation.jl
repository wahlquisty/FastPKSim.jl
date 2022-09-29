# Date: 220629
# Here, we assume that we are only looking for the concentration of the central compartment volume

using StaticArrays, SLEEFPirates

####

# update state for bolus
@inline function bolus(x, v)
    return x .+ v
end

# update state with input step
@inline function step(x, λ, Φdiag, u)
    Γdiag = @. 1 / λ * (Φdiag - 1)
    return @. Φdiag * x + Γdiag * u
end

# output computation
function gety(x, V1, R)
    return (1/V1*R'*x)[1]
end

# Diagonal of discrete time system matrix
@inline @fastmath function getΦdiag(λ, h)
    Φdiag = @. SLEEFPirates.exp(λ * h) # cannot differentiate through VectorizationBase.vexp
    # return @. exp(λ * h) 
    return Φdiag
end

# Initiate/update state
@inline function updatestate(x::AbstractVector{T}, h, λ, u=zero(T), v=zero(T)) where {T}
    Φdiag = getΦdiag(λ, h) # compute Φ
    x = bolus(x, v) # update state for bolus
    x = step(x, λ, Φdiag, u) # infusion affect next sample
    return x
end

# Update state and compute output
@inline function updatestateoutput(x::AbstractVector{T}, h, V1, λ, R, u=zero(T), v=zero(T)) where {T}
    Φdiag = getΦdiag(λ, h) # compute Φ
    x = bolus(x, v) # update state for bolus
    y = gety(x, V1, R) # compute output
    x = step(x, λ, Φdiag, u) # infusion affect next sample
    return x, y
end


"""
    pksim!(y, θ, u, v, hs, youts)
Fast simulation of the three compartment mammillary PK model.

The parameter vector θ has the following structure
```
θ = [k10, k12, k13, k21, k31, V1]
```
# Arguments:
- `y`: Preallocated output vector of size length(youts)
- `θ`: Parameter vector, see above.
- `u`: Infusion rate vector of size length(hs)
- `v`: Bolus dose vector of size length(hs)
- `hs`: Step size, should have the size of [diff(time) diff(time)[end]] where time is the matching time vector to u, v
- `youts`: Indices for output observations, corresponding to times in hs

Updates `y` with simulated outputs `x_1` at time instances `youts`.
"""
function pksim!(y, θ, u, v, hs, youts)
    λ, R = update(θ) # Setting up simulator
    j = 1 # counter to keep track of next free spot in y
    x = @SVector zeros(eltype(u), 3) # initial state
    for i in eachindex(u, hs, v)
        if i in youts # if we want to compute output
            x, yi = @inbounds updatestateoutput(x, hs[i], θ[6], λ, R, u[i], v[i]) # update state and compute output
            y[j] = yi
            j += 1
        else
            x = @inbounds updatestate(x, hs[i], λ, u[i], v[i]) # update state
        end
    end
    return y
end