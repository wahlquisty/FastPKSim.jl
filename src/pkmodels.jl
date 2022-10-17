# Mamillary three-compartment model functions
# Computes eigenvalues and R column vector of the three compartment mammillary model.

using StaticArrays, LinearAlgebra


struct PK{T,V,L,LI,S,O}
    θ::T
    V1inv::V
    λ::L
    λinv::LI
    R::S
    order::O
end

function PK(θ, order) # constructor
    λ = getλ(θ, order)
    R = getR(θ, λ, order)
    V1inv = 1 / θ[end]
    λinv = 1 ./ λ
    PK(θ, V1inv, λ, λinv, R, order)
end


function getλ(θ::AbstractVector{T}, order) where {T}
    if order == 1
        return -θ[1]

    elseif order == 2
        k10, k12, k21, _ = θ
        a1 = (k10 + k12 + k21) / 2
        a2 = sqrt(a1^2 - k10 * k21)
        return @SVector [-a1 - a2, -a1 + a2] # The eigenvalues of the continuous-time system matrix

    elseif order == 3

        k10, k12, k13, k21, k31, _ = θ
        b1 = k10 + k12 + k13 + k21 + k31
        b2 = k21 * (k10 + k13 + k31) + k31 * (k10 + k12)
        b3 = k10 * k21 * k31

        # Wengert list used to compute λ.
        a1 = b1 / 3
        a2 = a1^2
        a3 = a1 * a2
        a4 = b2 / 3
        a5 = a4 - a2
        a6 = (b1 * a4 - b3) / 2
        a7 = 2(a6 + sqrt(complex(a5^3 + (a3 - a6)^2)) - a3)^(1 / T(3))
        a8 = -real(a7)
        a9 = imag(a7)
        a10 = a9 * sqrt(T(3)) / 2
        a11 = a1 - a8 / 2

        return @SVector [-a1 - a8, -a10 - a11, a10 - a11] # The eigenvalues of the continuous-time system matrix
    end
end

# Computation of R, nominators of the first-order systems
@inline function getR(θ, λ, order) # samma för R
    if order == 1
        return one(1)
    elseif order == 2
        _, _, k21, _ = θ
        l1, l2 = λ
        a1 = l1 - l2
        a1inv = 1 / a1
        Qinv = @SMatrix [l1*a1inv -a1inv; -l2*a1inv a1inv]
        b = @SVector [1, k21] # First column of P, only interested in the first output, x1
        return Qinv * b
    elseif order == 3
        _, _, _, k21, k31, _ = θ
        l1, l2, l3 = λ
        a1 = l2 - l1
        a2 = l3 - l1
        a3 = l3 - l2
        d1 = a2 * a1
        d2 = -a3 * a1
        d3 = a2 * a3
        d1inv = 1 / d1
        d2inv = 1 / d2
        d3inv = 1 / d3
        Qinv = @SMatrix [(l1*d1inv)*l1 l1*d1inv d1inv; (l2*d2inv)*l2 l2*d2inv d2inv; (l3*d3inv)*l3 l3*d3inv d3inv]
        b = @SVector [1, k21 + k31, k21 * k31] # Quite often we would only be interested in the first column (first output). See paper for computations.
        return Qinv * b
    end
end