using Plots, SparseArrays, LinearAlgebra, Printf, LaTeXStrings

η(x) = .25*((1+tanh(20(x-.25)))*((1+tanh(20(-x+.75)))))
κ(x) = (sin(2π*x-π))^4

function CN(U,k,κ)
    m = length(U) - 1
    h = 1/(m+1)
    x = h:h:1 # include right end point
    A = Tridiagonal(fill(-1.0,m),fill(2.0,m+1),fill(-1.0,m)) |> sparse
    Aₖ = (1/h^2)*spdiagm(κ.(x))*A
    r = k/2
    U = (I+r*Aₖ)\((I-r*Aₖ)*U)
end

function RK2(U,k)
    f(u) = u*(1-u)*(u-1/2)
    U₀ = U + k*f.(U)/2
    U .+= k*f.(U₀)
end

h = 0.01
k = h

T = 1.
n = convert(Int64,ceil(T/k))
x = h:h:1 # include right end point
U = η.(x)
t = 0.
for i = 2:n+1
    global t += k
    global U = RK2(U,k/2)
    global U = CN(U,k,κ)
    global U = RK2(U,k/2)
    if t ≈ .01 || t ≈ .1 || t ≈ 1.
        p = plot(x, U, legend=:false)
        title!(@sprintf("Strang splitting solution at t = %1.2f",t))
        xlabel!(L"x")
        ylabel!(L"u(x,t)")
        display(p)
        savefig(p,@sprintf("p4at%1.2f.pdf",t))
    end
end

#repeat for cosine
κ(x) = (cos(2π*x-π))^4
U = η.(x)
t = 0.
for i = 2:n+1
    global t += k
    global U = RK2(U,k/2)
    global U = CN(U,k,κ)
    global U = RK2(U,k/2)
    if t ≈ .01 || t ≈ .1 || t ≈ 1.
        p = plot(x, U, legend=:false)
        title!(@sprintf("Strang splitting solution at t = %1.2f",t))
        xlabel!(L"x")
        ylabel!(L"u(x,t)")
        display(p)
        savefig(p,@sprintf("p4bt%1.2f.pdf",t))
    end
end
