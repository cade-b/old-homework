using Plots, Printf, LaTeXStrings, LinearAlgebra
λ, σ, μ = 10., 200., 100. #constants

import Base.isless #define complex ordering
function isless(a::ComplexF64,b::ComplexF64)
    return imag(a) < imag(b)
end

gfunc = u -> u*(1-u)^2
gprime = u -> 3u^2-4u+1 #derivative of g
f = u -> [u[4]-μ*u[3]+λ*gfunc(u[1]+u[2]), 
    μ*u[3]-u[4]+λ*gfunc(u[1]+u[2]), -σ*u[4], σ*u[3]]

η = [0.1, 0.1, 0., 0.1]

k, T = 0.0001, 20.

#define Jacobian
Df = u -> [λ*gprime(u[1]+u[2]) λ*gprime(u[1]+u[2]) -μ 1.;
    λ*gprime(u[1]+u[2]) λ*gprime(u[1]+u[2]) μ 1.;
    0. 0. 0. -σ;
    0. 0. σ 0.]

#Newton's method from demo
function Newton(x,g,Dg; tol = 1e-13, nmax = 100)
        for j = 1:nmax
            idk = Dg(x)
            dka = g(x)
            step = Dg(x)\g(x)
            x -= step
            if maximum(abs.(step)) < tol
                break
            end
            if j == nmax
                println("Newton's method did not terminate")
            end
        end
        x
end

#Trapezoid
g_tr = (u,Un,t,tn) -> u - Un - (k/2)*(f(u)+f(Un))
Dg_tr = u -> I - (k/2)*Df(u)
n = convert(Int64,ceil(T/k))
U = zeros(n+1,length(η)) # To save the solution values
U[1,:] = η
t = zeros(n+1)
t[1] = 0.
max_iter = 10
for i = 2:n+1
    t[i] = t[i-1] + k
    U[i,:] = Newton(U[i-1,:],u -> g_tr(u,U[i-1,:],t[i],t[i-1]), Dg_tr; 
        tol = k^3/10) 
end

p1 = plot(t, U,label = [L"u_1(t)" L"u_2(t)" L"u_3(t)" L"u_4(t)"])
xlabel!(L"t")
title!("\"True\" solution given by trapezoid")
display(p1)
savefig(p1,"trap.pdf")

λs = map( t -> eigvals(Df(U[convert(Int64,round(t/k+1)),:])), t)
λs = [ sort(i) for i in λs]
λ1 = [i[1] for i in λs]
λ2 = [i[2] for i in λs]
λ3 = [i[3] for i in λs]
λ4 = [i[4] for i in λs]

p2 = scatter(λ1 |> real, λ1 |> imag, markercolor = :blue, markerstrokewidth=0,
    label = L"\lambda_1")
scatter!(λ2 |> real, λ2 |> imag, markercolor = :green, markerstrokewidth=0, 
    label = L"\lambda_2")
scatter!(λ3 |> real, λ3 |> imag, markercolor = :red, markerstrokewidth=0, 
    label = L"\lambda_3")
scatter!(λ4 |> real, λ4 |> imag, markercolor = :yellow, markerstrokewidth=0, 
    label = L"\lambda_4")
xlabel!(L"\Re(\lambda_i)")
ylabel!(L"\Im(\lambda_i)")
title!("Eigenvalue distribution of Jacobian")
display(p2)
savefig(p2,"eigvals.pdf")

#Forward Euler
k, T = 0.001, 0.25
n = convert(Int64,ceil(T/k))
U_fe = zeros(n+1,length(η)) # To save the solution values
U_fe[1,:] = η
t = zeros(n+1) # To save times
t[1] = 0.
for i = 2:n+1
    U_fe[i,:] = U_fe[i-1,:] + k*f(U_fe[i-1,:])
    t[i] = t[i-1] + k
end
p3 = plot(t, U_fe,label = [L"u_1(t)" L"u_2(t)" L"u_3(t)" L"u_4(t)"])
xlabel!(L"t")
title!("Forward Euler solution")
display(p3)
savefig(p3,"fe.pdf")

#Leapfrog🐸
k, T = 0.001, 18.09
n = convert(Int64,ceil(T/k))
U_🐸 = zeros(n+1,length(η)) # To save the solution values
U_🐸[1,:] = η
t = zeros(n+1) # To save times
t[1] = 0.
#step of Forward Euler
U_🐸[2,:] = U_🐸[1,:] + k*f(U_🐸[1,:])
t[2] = k
for i = 3:n+1
    U_🐸[i,:] = U_🐸[i-2,:] + (2*k)*f(U_🐸[i-1,:]) #Leapfrog
    t[i] = t[i-1] + k
end
p4 = plot(t, U_🐸,label = [L"u_1(t)" L"u_2(t)" L"u_3(t)" L"u_4(t)"])
xlabel!(L"t")
title!("Leapfrog solution")
display(p4)
savefig(p4,"lf.pdf")