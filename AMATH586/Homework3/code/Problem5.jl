using Plots, LaTeXStrings

#build stability polynomial
function R_N(z,N)
    eval = 0
    for j = 0:N+1
        eval += z^j
    end
    eval
end

xrange = [-3,3]; yrange = [-2,2]
for N = [2 5 10 20 50]
    p = contourf(xrange[1]:0.01(1+rand()/10):xrange[2],
        yrange[1]:0.01(1+rand()/10):yrange[2],
        (x,y)-> sign(-(abs(R_N(x+1im*y,N)))+1),colorbar=false)
    xlabel!(L"\Re(z)")
    ylabel!(L"\Im(z)")
    title!(latexstring(@sprintf("\\mathrm{Stability~region~for}~N=%i",N)))
    display(p)
    savefig(p,@sprintf("N%i.pdf",N))
end

