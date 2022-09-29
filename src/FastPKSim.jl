module FastPKSim

export main

# Write your package code here.
function main()
    n = 10
    a = zeros(n)
    for i = 1:n
        a[i] = squareroot(i)
    end
    return a
end

function squareroot(x)
    return sqrt(x)
end

end