using Random, Plots, Distributions, LaTeXStrings, Printf

#part b
function simul(a, tmax, λ, μ, seed=123)
    Random.seed!(seed) #set seed

    t, n = 0., a
    values = [t a]
    while t<=tmax
        #time until next event
        dist1 = Exponential(1/(n*(λ+μ))) #Reciprocal of input sets exponential parameter
        T = rand(dist1)
        t += T

        #determine if event is birth or death
        dist2 = Bernoulli(λ/(λ+μ))
        result = 2rand(dist2)-1 #rescale to take values -1 and 1
        n += result
        values = [values; [t n]]

        #safegaurd to prevent dividing by 0
        if n == 0
            return values
            break
        end
    end
    values
end

#part c
#=λ, μ, a, t = 1., 1., 10, 100.
p1 = plot(title=latexstring(@sprintf(
    "\\mathrm{Birth-Death~Process~with}~\\lambda = %1.0f,~\\mu=%1.0f,a=%i,t=%i"
    ,λ, μ, a, t)),legend=false)
xlabel!(latexstring(@sprintf("\\mathrm{time}")))
ylabel!(latexstring(@sprintf("\\mathrm{Number~of~individuals}")))
for i=1:10
    results=simul(a, t, λ, μ, i)
    plot!(results[:,1],results[:,2])
end
display(p1)

#part d
λ, μ, a, t = 1.1, 1., 10, 50.
p2 = plot(title=latexstring(@sprintf(
    "\\mathrm{Birth-Death~Process~with}~\\lambda = %1.1f,~\\mu=%1.0f,a=%i,t=%i"
    ,λ, μ, a, t)),legend=false, yaxis=:log)
xlabel!(latexstring(@sprintf("\\mathrm{time}")))
ylabel!(latexstring(@sprintf("\\mathrm{Number~of~individuals}")))
for i=1:10
    results=simul(a, t, λ, μ, i)
    # remove zero value from log-scale plot
    if results[end,2]==0
        results = results[1:(end-1),:] 
    end
    plot!(results[:,1],results[:,2])
end
display(p2)=#

#part e
function onerun(t, a, λ, μ, seed=123)
    results = simul(a, t, λ, μ, seed)
    index = maximum(findall(results[:,1].<t))
    results[index,2]
end

function simulation(t,numruns, a, λ, μ, seed=123)
    total = 0
    for i=1:numruns
        total += onerun(t, a, λ, μ, seed)
    end
    total/numruns
end

#part f1
λ, μ, a = 1., 1., 10
tvec = 0:10:100
for (i,t) in enumerate(tvec)
    
