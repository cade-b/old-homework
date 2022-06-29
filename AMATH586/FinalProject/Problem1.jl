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
    α, ϵ = 0.5, 0.01

    g(v) = v*(α-v)*(v-1)-w
    f(v) = g(v)/ϵ
    fprime(v) = ForwardDiff.derivative(f,v) #get derivative with autodiff

    G₁ = (u,Un) -> u - Un - (k/4)*(f(Un)+f(u))
    dG₁ = u -> 1-(k/4)*fprime(u) #derivative of G₁
    U₀ = Newton(U,u -> G₁(u,U), dG₁)

    G₂ = (u,Ustar,Un) -> u-(1/3)*(4Ustar-Un+k*f(u))
    dG₂ = u -> 1-(k/3)*fprime(u)

    Newton(U,u -> G₂(u,U₀,U), dG₂)
end

w = 0.1
T = 0.1
k = 1e-5 #chosen to be small 
n = convert(Int64,ceil(T/k))
V = zeros(n+1)
V[1] = 0.3 #inital condition
for i = 2:n+1
    V[i] = TRBDF2(V[i-1], w, k)
end

p = plot(0:k:T,V, label=:false)
xlabel!(L"t")
ylabel!(L"v(t)")
title!(@sprintf("Computed solution with k=%0.0e",k))
savefig(p, "prob1.pdf")
display(p)


error = Inf
println("    t      |      error       | reduction ratio ")
for j = 4:10
    error_prev = error
    local k = 2.0^(-j)*0.1 #choose various values
    local n = convert(Int64,ceil(T/k))
    Vapprox = zeros(n+1)
    Vapprox[1] = 0.3 #inital condition
    for i = 2:n+1
        Vapprox[i] = TRBDF2(Vapprox[i-1], w, k)
    end
    global error = abs(V[end]-Vapprox[end])
    error_diff = error_prev/error
    println(@sprintf("%0.4e | %1.10e | %1.10e",k,error,error_diff))
end