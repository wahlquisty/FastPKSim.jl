# PK simulation example

using FastPKSim, Plots

θ = [1.f0,2.f0,3.f0,4.f0,5.f0,6.f0] # Parameter values, k10, k12, k13, k21, k31, V1

nu = 50 # length of input
u = [1.1f0*ones(Float32,10) ; 0.1f0*ones(Float32,nu-10)] # Infusion rates
v = [2.0f0*ones(Float32,10) ;0.0f0*ones(Float32,nu-10)  ] # Bolus doses

time = 0.f0:1.f0:nu
hs = diff(time)

youts = [2, 3, 6, 7, 9, 13, 16, 27, 46, 47] # observation times, matching time vector

y = zeros(length(youts)) # Create output vector
pksim!(y, θ, u, v, hs, youts) # Simulate model

plot(time[youts],y,xlabel="time",ylabel="y",label="")

# savefig("simulation_example.png")