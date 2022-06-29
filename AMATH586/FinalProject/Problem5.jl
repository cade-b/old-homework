using ForwardDiff, LinearAlgebra, Printf, Plots, LaTeXStrings

#Newton's method from demo
function Newton(x,g,Dg; tol = 1e-8, nmax = 1000)
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

#function created in problem 2
function TRBDF2_scalar(U,w,Iₐ,k)
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

#use map to vectorize scalar TR-BDF-2 function
function TRBDF2(V,W,k)
    map((v,w,iₐ) -> TRBDF2_scalar(v,w,iₐ,k),V,W,Iₐvec) 
end

function RK2(V,W,k)
    β, γ = 1., 1.

    f(W,V) = β*V-γ*W
    W₀ = W + (k/2)*f(W,V)
    W += k*f(W₀,V)
end

function TRBDF2_heat(V)
    V₀ = (I-(k₂/4)*B)\((I+(k₂/4)*B)*V)
    (I-(k₂/3)*B)\((1/3)*(4V₀-V))
end

h = 0.02
k = h/10
global k₂ = k/2 #get half timestep
a, b, κ = 0., 6., 0.2
x = a:h:b
m = length(x)-2
T = 4.
n = convert(Int64,ceil(T/k))

A = Tridiagonal(fill(1.0,m+1),fill(-2.0,m+2),fill(1.0,m+1))
A[1,2] = 2.
A[end,end-1] = 2.
global B = (κ/h^2)*A #add scaling

V = zeros(m+2)
W = zeros(m+2)

anim = Animation()
plot(x,V, yaxis = [-0.1,1.1], label=L"v(x,t)") 
plot!(x,W, yaxis = [-0.1,1.1], label=L"w(x,t)")
xlabel!(L"x")
title!(latexstring("t=0"))
frame(anim)
savefig("prob5_t=0.pdf")

global Iₐvec = @. 0.8exp(-5x^2)

for i=2:n+1
    global V = TRBDF2_heat(V)
    global W = RK2(V,W,k/2)
    V = TRBDF2(V,W,k)
    W = RK2(V,W,k/2)
    V = TRBDF2_heat(V)
    if mod(i,10)==0
        plot(x,V, yaxis = [-0.1,1.1], label=L"v(x,t)") 
        plot!(x,W, yaxis = [-0.1,1.1], label=L"w(x,t)")
        xlabel!(L"x")
        title!(latexstring(@sprintf("t=%1.2f",i*k)))
        frame(anim)
        if mod(i,500)==0
            savefig(@sprintf("prob5_t=%1.0f.pdf",i*k))
        end
    end
end
gif(anim,"prob5.gif")