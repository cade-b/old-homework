using Plots, Distributions, Random, LaTeXStrings, Printf

u, T, s = 3.4e-5, 3/365, 0.004
keepExtinct = false
Tₑ = 25. #final time in years
Random.seed!(1234)

function simul(u, T, s, Tₑ, keepExtinct)
    n = convert(Int64,ceil(Tₑ/T))

    x₁ = zeros(Int64,n+1)
    x₂ = zeros(Int64,n+1)
    t = zeros(n+1)
    x₁[1], x₂[1], t[1] = 1, 0, 0.

    d₁ = (1-s)/2
    d₂ = (1-s)^2/2
    b₁ = 1-d₁
    b₂ = 1-d₂
    probs₁ = [b₁*(1-u), d₁, b₁*u]
    probs₂ = [b₂, d₂, 0]

    for i = 2:n+1
        dist₁ = Multinomial(x₁[i-1], probs₁)
        dist₂ = Multinomial(x₂[i-1], probs₂)
        result₁ = rand(dist₁)
        result₂ = rand(dist₂)
        x₁[i] = x₁[i-1]+result₁[1]-result₁[2]
        x₂[i] = x₂[i-1]+result₂[1]-result₂[2]+result₁[3]
        t[i] = t[i-1]+T

        if (x₁[i]==0 && x₂[i]==0) && keepExtinct == false
            return simul(u, T, s, Tₑ, keepExtinct)
        end
    end
    return x₁, x₂, t
end

trial1₁, trial1₂, t = simul(u, T, s, Tₑ, keepExtinct)
ind₁ = trial1₁.>0
ind₂ = trial1₂.>0
p1 = plot(t[ind₁],trial1₁[ind₁],label=L"N_1(t)", yaxis=:log)
plot!(t[ind₂],trial1₂[ind₂],label=L"N_2(t)")
xlabel!("time (years)")
title!("Simulation 1")
display(p1)
savefig(p1,"sim1.pdf")

trial2₁, trial2₂, t = simul(u, T, s, Tₑ, keepExtinct)
ind₁ = trial2₁.>0
ind₂ = trial2₂.>0
p2 = plot(t[ind₁],trial2₁[ind₁],label=L"N_1(t)", yaxis=:log)
plot!(t[ind₂],trial2₂[ind₂],label=L"N_2(t)")
xlabel!("time (years)")
title!("Simulation 2")
display(p2)
savefig(p2,"sim2.pdf")

trial3₁, trial3₂, t = simul(u, T, s, Tₑ, keepExtinct)
ind₁ = trial3₁.>0
ind₂ = trial3₂.>0
p3 = plot(t[ind₁],trial3₁[ind₁],label=L"N_1(t)", yaxis=:log)
plot!(t[ind₂],trial3₂[ind₂],label=L"N_2(t)")
xlabel!("time (years)")
title!("Simulation 3")
display(p3)
savefig(p3,"sim3.pdf")

#part b
function compute_avg(u, T, s, Tₑ, numruns)
    n = convert(Int64,ceil(Tₑ/T))
    avg₁ = zeros(n+1)
    avg₂ = zeros(n+1)
    for i = 1:numruns
        x₁, x₂, t = simul(u, T, s, Tₑ, true)
        avg₁ .+= x₁./numruns
        avg₂ .+= x₂./numruns
    end
    return avg₁, avg₂, t
end

numruns = 100000 #chosen number of runs
avg₁, avg₂, t = compute_avg(u, T, s, Tₑ, numruns)
ind₁ = avg₁.>0
ind₂ = avg₂.>0
p4 = plot(t[ind₁],avg₁[ind₁],label=L"N_1(t)", yaxis=:log)
plot!(t[ind₂],avg₂[ind₂],label=L"N_2(t)")
xlabel!("time (years)")
title!(@sprintf("Average over %i runs",numruns))
display(p4)
savefig(p4,"avg.pdf")

#part c
d₁ = (1-s)/2
d₂ = (1-s)^2/2
b₁ = 1-d₁
b₂ = 1-d₂
xₜ = t -> (1+b₁*(1-u)-d₁)^(round(t/T))
yₜ = t -> (b₁*u*((1+b₂-d₂)^(round(t/T))-(1+b₁*(1-u)-d₁)^(round(t/T))))/((b₂-d₂)-(b₁*(1-u)-d₁))

ind₁ = xₜ.(t).>0
ind₂ = yₜ.(t).>0
p5 = plot(t[ind₁], xₜ.(t[ind₁]),label=L"N_1(t)", yaxis=:log)
plot!(t[ind₂],yₜ.(t[ind₂]),label=L"N_2(t)")
xlabel!("time (years)")
title!("Average number of cells computed analytically")
display(p5)
savefig(p5,"true_avg.pdf")
