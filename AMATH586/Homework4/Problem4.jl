using LinearAlgebra, Plots, SparseArrays, LaTeXStrings, Random, UncertainData

Random.seed!(123)

function prand(m)
    p = x -> -(2.0/3)*x .+4.0/3 .+ .5sin.(2*pi*x)
    B = 1.7
    out = fill(0., m)
    for j = 1:m
        u = 10.
        y = 0.
        while u >= p(y)/ B
            y = rand()
            u = rand()
        end
        out[j] = y
    end
    out
end

h, s = 0.0001, 2.
k = 10*h
T = 0.1 #final time

x = 0:h:1
m = length(x)-2
N = m

#build Yᵢ
X = prand(N)
bins = bin(x[2:end], X, ones(length(X)))
Y = sum.(bins)

x = x[2:end-1] #remove BC

#build initial condition
u₀ = Y./(h*N)

#code from problem 2
B = spdiagm(-1 => ones(m-1), 0 => -2*ones(m), 1 => ones(m-1))
B[1,1] += s/(1+s)
B[1,end] += s/(1+s)
B[end,1] += 1/(1+s)
B[end,end] += 1/(1+s)

A = I - (k/(2h^2))*B

n = convert(Int64,ceil(T/k))
U = zeros(m,n+1)
U[:,1] = u₀ 
t = zeros(n+1)
t[1] = 0
for i = 2:n+1
    t[i] = t[i-1] + k
    U[:,i] = A\U[:,i-1]
end

ρ = x -> -2x/3+4/3+sin(2π*x)/2
p1 = plot(x, ρ.(x), label=L"\rho(x)", linewidth=3)
xlabel!("x")
title!("Density approximations versus true density")

ind₁ = t.≈ 0.001
plot!(x,U[:,ind₁], label=L"t=0.001")

ind₂ = t.≈ 0.01
plot!(x,U[:,ind₂], label=L"t=0.01")

ind₃ = t.≈ 0.1
plot!(x,U[:,ind₃], label=L"t=0.1")
savefig(p1, "p4.pdf")
display(p1)