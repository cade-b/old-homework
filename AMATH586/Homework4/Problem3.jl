using Random, Plots, LaTeXStrings, Printf

Random.seed!(123)

n = 10000
X = randn(n)
x = -3:0.001:3

ρ = x -> exp(-x^2/2)/√(2π)

function density_approx(x, t, X)
    N = length(X)
    eval = zeros(length(x))
    for j = 1:N
       eval += @. exp(-(x-X[j])^2/(2t))/√(2π*t)
    end
    eval ./= N
end

p1 = plot(x, ρ.(x), label=L"\rho(x)", linewidth=3)
xlabel!("x")
title!("Density approximations versus (old) true density")

for t = [0.001 0.01 0.1 1. 10.]
    plot!(x, density_approx(x, t, X),label=@sprintf("t=%2.3f",t))
end
display(p1)
savefig(p1,"density_approx.pdf")

# repeating with new density
function prand(m)
    p = x -> -(2.0/3)*x .+4.0/3 .+ .5sin.(2*pi*x)
    B = 1.7
    out = fill(0., m)
    for j = 1:m
        u = 10.
        y = 0.
        while u >= p(y)/ B
            y = rand()
            u = rand()
        end
        out[j] = y
    end
    out
end


x = 0:0.001:1 #new domain where density should be nonzero
ρ = x -> -2x/3+4/3+sin(2π*x)/2
X = prand(n)
p2 = plot(x, ρ.(x), label=L"\rho(x)", linewidth=3)
xlabel!("x")
title!("Density approximations versus (new) true density")

for t = [0.001 0.01 0.1 1. 10.]
    plot!(x, density_approx(x, t, X),label=@sprintf("t=%2.3f",t))
end
display(p2)
savefig(p2,"density_approx_new.pdf")