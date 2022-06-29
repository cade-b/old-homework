using LinearAlgebra, Printf, Plots, LaTeXStrings


function TRBDF2_heat(V)
    V₀ = (I-(k/4)*B)\((I+(k/4)*B)*V)
    (I-(k/3)*B)\((1/3)*(4V₀-V))
end

h = 0.01
global k = h #global to avoid passing
a, b, κ = 0., 4., 1.
x = a:h:b
m = length(x)-2
T = 6.
n = convert(Int64,ceil(T/k))

A = Tridiagonal(fill(1.0,m+1),fill(-2.0,m+2),fill(1.0,m+1))
A[1,2] = 2.
A[end,end-1] = 2.
global B = (κ/h^2)*A #add scaling

V = zeros(m+2,n+1)
#implement initial condition
ind = x.<(a+b)/2
V[ind,1] .= 1

anim = Animation()
p = plot(x,V[:,1], yaxis = [-0.1,1.1], label=:false)
xlabel!(L"x")
ylabel!(L"v(x,t)")
title!(latexstring("t=0"))
frame(anim)
savefig("prob3_t=0.pdf")

for i=2:n+1
    V[:,i] = TRBDF2_heat(V[:,i-1])
    plot(x,V[:,i], yaxis = [-0.1,1.1], label=:false)
    if mod(i,2)==0
        xlabel!(L"x")
        ylabel!(L"v(x,t)")
        title!(latexstring(@sprintf("t=%1.2f",i*k)))
        frame(anim)
        if i*k ≈ round(i*k)
            savefig(@sprintf("prob3_t=%1.0f.pdf",i*k))
        end
    end
end
gif(anim,"prob3.gif")

avgV = h*sum(V,dims=1)
p = plot(0:k:T,-avgV'.+2, label=:false)
xlabel!(L"t")
title!("Difference between true and computed integrals")
savefig(p,"prob3_discrep.pdf")
display(p)

#try smoother initial condition
V2 = zeros(m+2,n+1)
sigmoid(x) = 1-1/(1+exp(-10(x-2)))
V2[:,1] = sigmoid.(x)
for i=2:n+1
    V2[:,i] = TRBDF2_heat(V2[:,i-1])
end

avgV2 = h*sum(V2,dims=1)
p1 = plot(0:k:T,abs.(avgV2'.-avgV2[1]), label=:false)
xlabel!(L"t")
title!("Heat loss for sigmoid")
savefig(p1,"prob3_sigmoid.pdf")
display(p1)