using Plots, SparseArrays, LinearAlgebra, Printf, ForwardDiff, LaTeXStrings

f(x) = exp(-20(x-1/2)^2)
fprime(x) = ForwardDiff.derivative(f,x) #get initial x derivative with autodiff
g(x) = sin(4π*x)

c = 1.0
h = 0.01
k = h #choose k=h for stability

m = convert(Int64,1/h)-1
T = 3.
n = convert(Int64,ceil(T/k))

A = Tridiagonal(fill(-1.0,m),fill(2.0,m+1),fill(-1.0,m)) |> sparse
A[1,end] = -1
A[end,1] = -1

S = Tridiagonal(fill(-1.0,m),fill(0.0,m+1),fill(1.0,m)) |> sparse
S[1,end] = -1
S[end,1] = 1

A1 = I - (-c*k)/2h*S - (k^2*(-c)^2)/2h^2*A; # Lax-Wendroff matrix for 1st step
A2 = I - (c*k)/2h*S - (k^2*c^2)/2h^2*A; # Lax-Wendroff matrix for 2nd step

x = h:h:1 # include right end point
u = zeros(m+1, n+1)
w = zeros(m+1, n+1)
u[:,1] = f.(x)
w[:,1] = g.(x) -c*fprime.(x)
t = 0.
p = plot(x, u[:,1])
title!("Lax-Wendroff solution at t = 0", legend=:false)
        xlabel!(L"x")
        ylabel!(L"u(x,t)")
        display(p)
        savefig(p,"p3t0.pdf")

for i = 2:n+1
    global t += k
    u[:,i] = A1*u[:,i-1] + k*w[:,i-1]
    w[:,i] = A2*w[:,i-1]
    if t ≈ round(t)
        p = plot(x, u[:,i])
        title!(@sprintf("Lax-Wendroff solution at t = %i",t), legend=:false)
        xlabel!(L"x")
        ylabel!(L"u(x,t)")
        display(p)
        savefig(p,@sprintf("p3t%i.pdf",t))
    end
end


