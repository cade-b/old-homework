using ForwardDiff, LinearAlgebra, Printf, Plots, LaTeXStrings

#Newton's method from demo
function Newton(x,g,Dg; tol = 1e-13, nmax = 1000)
    for j = 1:nmax
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

function TRBDF2(U,w,k)
    α, ϵ = 0.3, 0.001

    g(v) = v*(α-v)*(v-1)-w+Iₐ
    f(v) = g(v)/ϵ
    fprime(v) = ForwardDiff.derivative(f,v) #get derivative with autodiff

    G₁ = (u,Un) -> u - Un - (k/4)*(f(Un)+f(u))
    dG₁ = u -> 1-(k/4)*fprime(u)
    U₀ = Newton(U,u -> G₁(u,U), dG₁)

    G₂ = (u,Ustar,Un) -> u-(1/3)*(4Ustar-Un+k*f(u))
    dG₂ = u -> 1-(k/3)*fprime(u)

    Newton(U,u -> G₂(u,U₀,U), dG₂)
end

function RK2(V,W,k)
    β, γ = 1., 1.

    f(W,V) = β*V-γ*W
    W₀ = W + k*f(W,V)/2
    W += k*f(W₀,V)
end

T = 1.
k = 1e-5
n = convert(Int64,ceil(T/k))
global Iₐ = 0.
V = zeros(n+1)
W = zeros(n+1)
V[1] = 0.31 #inital condition
W[1] = 0.
for i = 2:n+1
    W[i] = RK2(V[i-1],W[i-1],k/2)
    V[i] = TRBDF2(V[i-1], W[i], k)
    W[i] = RK2(V[i],W[i],k/2)
end

p1 = plot(0:k:T, V, label=L"v(t)")
plot!(0:k:T, W, label=L"w(t)")
xlabel!(L"t")
title!(L"\mathrm{Numerical~solution~with}~v_0=0.31")
savefig(p1,"prob2a.pdf")
display(p1)

T = .25
n = convert(Int64,ceil(T/k))
V = zeros(n+1)
W = zeros(n+1)
V[1] = 0.29 #new inital condition
W[1] = 0.
for i = 2:n+1
    W[i] = RK2(V[i-1],W[i-1],k/2)
    V[i] = TRBDF2(V[i-1], W[i], k)
    W[i] = RK2(V[i],W[i],k/2)
end

p2 = plot(0:k:T, V, label=L"v(t)")
plot!(0:k:T, W, label=L"w(t)")
xlabel!(L"t")
title!(L"\mathrm{Numerical~solution~with}~v_0=0.29")
savefig(p2,"prob2b.pdf")
display(p2)

error1, error2 = Inf, Inf
println("    t      |     error v      |     error w      
    |reduction ratio v |reduction ratio w")
for j = 4:10
    error1_prev = error1
    error2_prev = error2
    local k = 2.0^(-j)*0.1 #small 
    local n = convert(Int64,ceil(T/k))
    Va = zeros(n+1)
    Wa = zeros(n+1)
    Va[1] = 0.29 #inital condition
    Wa[1] = 0.
    for i = 2:n+1
        Wa[i] = RK2(Va[i-1],Wa[i-1],k/2)
        Va[i] = TRBDF2(Va[i-1], Wa[i], k)
        Wa[i] = RK2(Va[i],Wa[i],k/2)
    end
    global error1 = abs(V[end]-Va[end])
    global error2 = abs(W[end]-Wa[end])
    error1_diff = error1_prev/error1
    error2_diff = error2_prev/error2
    println(@sprintf("%0.4e | %1.10e | %1.10e | %1.10e | %1.10e
        ",k,error1, error2, error1_diff, error2_diff))
end

T = 10.
k = 0.001
n = convert(Int64,ceil(T/k))
global Iₐ = 0.2 #make global to avoid passing into TRBDF2
V = zeros(n+1)
W = zeros(n+1)
V[1] = 0. #inital condition
W[1] = 0.
for i = 2:n+1
    W[i] = RK2(V[i-1],W[i-1],k/2)
    V[i] = TRBDF2(V[i-1], W[i], k)
    W[i] = RK2(V[i],W[i],k/2)
end

p3 = plot(0:k:T, V, label=L"v(t)")
plot!(0:k:T, W, label=L"w(t)")
xlabel!(L"t")
title!(L"\mathrm{Numerical~solution~with}~I_a=0.2")
savefig(p3,"prob2c.pdf")
display(p3)