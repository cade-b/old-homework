using Plots, LaTeXStrings

R = z -> (12+5z)/((4-z)*(3-z))
xrange = [-5,15]; yrange = [-10,10]
p = contourf(xrange[1]:0.01(1+rand()/10):xrange[2],
    yrange[1]:0.01(1+rand()/10):yrange[2],(x,y)-> sign(-(abs(R(x+1im*y)))+1),
    colorbar=false)
xlabel!(L"\Re(z)")
ylabel!(L"\Im(z)")
title!("Stability region for TR-BDF2")
display(p)
savefig(p,"trbdf.pdf")
