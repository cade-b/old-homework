clear

d = domain(0,1);
x = chebfun('x',d);
utrue = x*(1-x);
f = 2*(3*x^2-x+1);
L = chebop(@(x,u) -diff((1+x^2)*diff(u)),d,0,0);
u = L\f;
err_2 = norm(u-utrue,2)
err_inf = norm(u-utrue,'inf')

