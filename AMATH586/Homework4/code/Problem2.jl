using LinearAlgebra, Plots, Printf, SparseArrays, LaTeXStrings

η = x -> exp(-20*(x-1/2)^2)

h, k = 0.001, 0.001
s = 2.
T = 0.1 #final time

x = 0:h:1
x = x[2:end-1] #remove BC
m = length(x)
B = spdiagm(-1 => ones(m-1), 0 => -2*ones(m), 1 => ones(m-1))
B[1,1] += s/(1+s)
B[1,end] += s/(1+s)
B[end,1] += 1/(1+s)
B[end,end] += 1/(1+s)

A = I - (k/(2h^2))*B

u₀ = η.(x)
n = convert(Int64,ceil(T/k))
U = zeros(m,n+1)
U[:,1] = u₀ 
t = zeros(n+1)
t[1] = 0
for i = 2:n+1
    t[i] = t[i-1] + k
    U[:,i] = A\U[:,i-1]
end

ind₁ = t.≈ 0.001
p1 = plot(x,U[:,ind₁], label=false)
xlabel!("x")
ylabel!("u(x,t)")
title!("Approximate solution at time 0.001")
savefig(p1, "p2_1.pdf")
display(p1)

ind₂ = t.≈ 0.01
p2 = plot(x,U[:,ind₂], label=false)
xlabel!("x")
ylabel!("u(x,t)")
title!("Approximate solution at time 0.01")
savefig(p2, "p2_2.pdf")
display(p2)

ind₃ = t.≈ 0.1
p3 = plot(x,U[:,ind₃], label=false)
xlabel!("x")
ylabel!("u(x,t)")
title!("Approximate solution at time 0.1")
savefig(p3, "p2_3.pdf")
display(p3)