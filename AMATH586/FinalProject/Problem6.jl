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
function TRBDF2_scalar(U,w,k)
    α, ϵ = 0.1, 0.01

    g(v) = v*(α-v)*(v-1)-w
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
    map((v,w) -> TRBDF2_scalar(v,w,k),V,W) 
end

function RK2(V,W,k)
    β, γ = 0.5, 1.

    f(W,V) = β*V-γ*W
    W₀ = W + (k/2)*f(W,V)
    W += k*f(W₀,V)
end

function TRBDF2_heat(V)
    V₀ = (I-(k₂/4)*B)\((I+(k₂/4)*B)*V)
    (I-(k₂/3)*B)\((1/3)*(4V₀-V))
end

h = 0.05
k = h/10
global k₂ = k/2 #get half timestep
a, b, κ = 0., 12., 1.
T = 0.9
n = convert(Int64,ceil(T/k))

x = a:h:b |> Array
y = x
X = repeat(reshape(x, 1, :), length(y), 1)
Y = repeat(reverse(y), 1, length(x));
Dₙ = (x,y) -> x < 2 ? 1.0 : 0.0
V = map(Dₙ,X,Y)
W = 0*V

m = length(x)-2
A = Tridiagonal(fill(1.0,m+1),fill(-2.0,m+2),fill(1.0,m+1))
A[1,2] = 2.
A[end,end-1] = 2.
global B = (κ/h^2)*A #add scaling

anim = Animation()
contourf(x,y,reverse(V,dims=1), clim=(-.5,1.5))
xlabel!(L"x")
ylabel!(L"y")
title!(latexstring("t=0"))
frame(anim)
savefig("prob6_t=0.pdf")

for i = 2:n+1
    global V = TRBDF2_heat(V)
    V = TRBDF2_heat(V')' |> Array
    global W = RK2(V,W,k/2)
    V = TRBDF2(V,W,k)
    W = RK2(V,W,k/2)
    V = TRBDF2_heat(V)
    V = TRBDF2_heat(V')' |> Array
    if mod(i,5)==0
        contourf(x,y,reverse(V,dims=1),clim=(-.5,1.5))
        xlabel!(L"x")
        ylabel!(L"y")
        title!(latexstring(@sprintf("t=%1.2f",i*k)))
        frame(anim)
        if mod(i,500)==0
            savefig(@sprintf("prob6_t=%1.0f.pdf",i*k))
        end
    end
end
savefig("prob6_t=0.9.pdf")

Dₙ = (x,y) -> y <= 6 
V = map(Dₙ,X,Y).*V #zero out necessary entries
T = 6.0-0.9
n = convert(Int64,ceil(T/k))
for i = 2:n+1
    global V = TRBDF2_heat(V)
    V = TRBDF2_heat(V')' |> Array
    global W = RK2(V,W,k/2)
    V = TRBDF2(V,W,k)
    W = RK2(V,W,k/2)
    V = TRBDF2_heat(V)
    V = TRBDF2_heat(V')' |> Array
    if mod(i,5)==0
        #need to flip matrix when plotting to account for ordering convention
        contourf(x,y,reverse(V,dims=1),clim=(-.5,1.5)) 
        xlabel!(L"x")
        ylabel!(L"y")
        title!(latexstring(@sprintf("t=%1.2f",i*k+0.9)))
        frame(anim)
        if i*k+0.9 ≈ round(i*k+0.9)
            savefig(@sprintf("prob6_t=%1.0f.pdf",i*k+0.9))
        end
    end
end

gif(anim,"prob6.gif")