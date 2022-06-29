using Plots, LaTeXStrings, Elliptic.Jacobi, Printf

β₁ = 0.
β₂ = 1.
β₃ = 10.
c = (β₁ + β₂ + β₃)/3
v = t -> β₂ + (β₃ - β₂)*cn(sqrt((β₃-β₁)/12)*t, (β₃-β₂)/(β₃-β₁) )^2

f = u -> [u[2], u[3], u[2]*(c - u[1])]
Df = u -> [0. 1. 0.; 0. 0. 1.; -u[2] c-u[1] 0.]
u₀ = [β₃,0.,-1.0/6*(β₃-β₁)*(β₃-β₂)]

h = 0.0001
[v(0), (v(h)-v(-h))/(2h), (v(h)-2v(0)+v(-h))/(h^2)]

T = 10.# Final time.

k = .02 #stepsize
p = 7
data = zeros(p)
ks = zeros(p)
for i = 1:p
    global k = k/2
    n = convert(Int64,ceil(T/k))
    println("Number of time steps = ", n)
    U = zeros(3,n+1) # To save the solution values
    U[:,1] = u₀
    global t = zeros(n+1,1)
    t[1] = 0.
    max_iter = 10
    # implementation of Heun's method
    for i = 2:n+1
        t[i] = t[i-1] + k
        Y1 = U[:,i-1]
        f1 = f(Y1)    
        Y2 = U[:,i-1] + (k/3)*f1
        f2 = f(Y2)
        Y3 = U[:,i-1] + (2k/3)*f2
        f3 = f(Y3)
        U[:,i] = U[:,i-1] + (k/4)*(f1+3*f3)
    end
    data[i] = abs(U[1,end] - v(t[end]))
    ks[i] = k
end
p1 = plot(ks,data,lw=2,ms=5,marker=:d, minorgrid = true, xaxis=:log, yaxis= :log,
    label=L"\mathrm{3rd~order~RK}", legend = :bottomright)
plot!(ks,ks.^3,lw=2,ms=5,marker=:d, minorgrid = true,
     label=L"\mathrm{O(k^3)~reference~line}")
xlabel!(latexstring("k"))
ylabel!(L"\mathrm{Error}")
title!("Convergence of RK3")
savefig(p1, "p1.pdf")
display(p1)

#print table
@printf("%s        | %s      | %s   \n","k","error","error reduction")
@printf("%f | %0.4e |       \n",ks[1],data[1])
for j=2:7
    @printf("%f | %0.4e |      %0.4f  \n",ks[j],data[j],data[j-1]/data[j])
end