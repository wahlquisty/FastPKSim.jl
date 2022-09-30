# Mamillary three-compartment model functions
# Computes eigenvalues and R column vector of the three compartment mammillary model.

using StaticArrays, LinearAlgebra

# Initiate/update parameters
function update(θ)
    λ = getλ(θ)
    λinv = 1 ./ λ
    R = getR(θ, λ)
    return λ, λinv, R
end

# Compute eigenvalues λ for 3 compartment mammillary model
@inline @fastmath function getλ(θ::AbstractVector{T}) where {T}
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

# Computation of R, nominators of the first-order systems
@inline function getR(θ, λ)
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