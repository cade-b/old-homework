using Plots, LaTeXStrings, LinearAlgebra, SparseArrays, Arpack, Printf

# 1st case
N = 6
T = spdiagm(-1=>0.5*ones(N-1),1=>0.5*ones(N-1))
initeig = eigs(T)

t = 100. # Final time.
k = 0.01 # Step size

# function to compute ST-TS
function rhsf(T)
    S = spdiagm(-1=>-diag(T,1),1=>diag(T,1))
    S*T-T*S
end

# implement RK4
n = convert(Int64,t/k) # Number of time steps, converted to Int64
tvec = k*(0:n)
A = zeros(N-1,n+1)
B = zeros(N,n+1)
A[:,1] = diag(T,1)
B[:,1] = diag(T)
for i = 2:n+1
    local Y1 = T
    f1 = rhsf(Y1)    
    Y2 = T + (k/2)*f1
    f2 = rhsf(Y2)
    Y3 = T + (k/2)*f2
    f3 = rhsf(Y3)
    Y4 = T + k*f3
    f4 = rhsf(Y4)
    global T += (k/6)*(f1+2*f2+2*f3+f4)
    A[:,i] = diag(T,1)
    B[:,i] = diag(T)
end

p1 = plot(
    title=latexstring("\\mathrm{RK4~Solution~for~b_j(t)~in~original~system}"))
xlabel!(latexstring("t"))
for j = 1:N
    plot!(tvec,B[j,:],label=latexstring(@sprintf("b_%i(t)",j)))
end
savefig(p1, "p7i.pdf")
display(p1)

p1a = plot(
    title=latexstring("\\mathrm{RK4~Solution~for~a_j(t)~in~original~system}"))
xlabel!(latexstring("t"))
for j = 1:N-1
    plot!(tvec,A[j,:],label=latexstring(@sprintf("a_%i(t)",j)))
end
savefig(p1a, "p7ii.pdf")
display(p1a)

# 2nd case
N = 12
T = spdiagm(0=>-2*ones(N),-1=>ones(N-1),1=>ones(N-1))
initeig2 = eigs(T)

# implement RK4
n = convert(Int64,t/k) # Number of time steps, converted to Int64
tvec = k*(0:n)
A2 = zeros(N-1,n+1)
B2 = zeros(N,n+1)
A2[:,1] = diag(T,1)
B2[:,1] = diag(T)
for i = 2:n+1
    local Y1 = T
    f1 = rhsf(Y1)    
    Y2 = T + (k/2)*f1
    f2 = rhsf(Y2)
    Y3 = T + (k/2)*f2
    f3 = rhsf(Y3)
    Y4 = T + k*f3
    f4 = rhsf(Y4)
    global T += (k/6)*(f1+2*f2+2*f3+f4)
    A2[:,i] = diag(T,1)
    B2[:,i] = diag(T)
end

p2 = plot(
    title=latexstring("\\mathrm{RK4~Solution~for~b_j(t)~in~new~system}"))
xlabel!(latexstring("t"))
for j = 1:N
    plot!(tvec,B2[j,:],label=latexstring(@sprintf("b_{%i}(t)",j)))
end
savefig(p2, "p7iii.pdf")
display(p2)

p2a = plot(
    title=latexstring("\\mathrm{RK4~Solution~for~a_j(t)~in~new~system}"))
xlabel!(latexstring("t"))
for j = 1:N-1
    plot!(tvec,A2[j,:],label=latexstring(@sprintf("a_{%i}(t)",j)))
end
savefig(p2a, "p7iv.pdf")
display(p2a)