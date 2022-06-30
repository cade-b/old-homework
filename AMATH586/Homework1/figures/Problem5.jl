using Plots, LaTeXStrings, Printf

α,β,γ,δ = 1.0,1.0,1.0,1.0 #parameters
T = 50. #final time
k = 0.001 #step size

f = u-> [α*u[1]-β*u[1]*u[2], δ*u[1]*u[2]-γ*u[2]] #function f such that u'(t)=f(u)
u0 = [5., 0.8] #initial condition

# Forward Euler
n = convert(Int64,T/k)# Number of time steps, converted to Int64
Uf = zeros(2,n+1) # To save the solution values
Uf[:,1] = u0
for i = 2:n+1
    Uf[:,i] = Uf[:,i-1] + k*f(Uf[:,i-1])
end
tvec = k*(0:n) #time vector

p1 = (plot(Uf[1,:],Uf[2,:], title = latexstring(@sprintf(
    "\\mathrm{Forward~Euler~with}~k = %1.3f",k)),legend=false))
xlabel!(latexstring("u_1(t)"))
ylabel!(latexstring("u_2(t)"))
savefig(p1, "fe1.pdf")
display(p1)
p2 = plot(tvec,Uf[1,:], title = latexstring(@sprintf(
    "\\mathrm{Forward~Euler~with}~k = %1.3f",k)),label=latexstring("u_1(t)"))
plot!(tvec,Uf[2,:],label=latexstring("u_2(t)"))
xlabel!(latexstring("t"))
savefig(p2, "fe2.pdf")
display(p2)

#Newton's method from sample notebook
function Newton!(x,g,Dg; tol = 1e-13, nmax = 100)
    for j = 1:nmax
        step = Dg(x)\g(x)
        x[1:end] -= step
        if maximum(abs.(step)) < tol
            break
        end
        if j == nmax
            println("Newton's method did not terminate")
        end
    end
    x
end

g = (U,Un) -> U - Un - k*f(U)
Dg = (U) -> [1.0+k*(α-β*U[2]) -k*β*U[1]; 
    k*δ*U[2] 1.0+k*(δ*U[1]-γ)]

# Backward Euler
n = convert(Int64,T/k) # Number of time steps, converted to Int64
Ub = zeros(2,n+1) # To save the solution values
Ub[:,1] = u0
for i = 2:n+1
    Unew = Ub[:,i-1] |> copy
    Newton!(Unew,u -> g(u,Ub[:,i-1]), Dg)    
    Ub[:,i] = Unew
end

p3 = plot(Ub[1,:],Ub[2,:], title = latexstring(@sprintf(
    "\\mathrm{Backward~Euler~with}~k = %1.3f",k)),legend=false)
xlabel!(latexstring("u_1(t)"))
ylabel!(latexstring("u_2(t)"))
savefig(p3, "be1.pdf")
display(p3)
p4 = plot(tvec,Ub[1,:], title = latexstring(@sprintf(
    "\\mathrm{Backward~Euler~with}~k = %1.3f",k)),label=latexstring("u_1(t)"))
plot!(tvec,Ub[2,:],label=latexstring("u_2(t)"))
xlabel!(latexstring("t"))
savefig(p4, "be2.pdf")
display(p4)