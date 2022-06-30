clear
ufunc = @(x) x.*(1-x);
f = @(x) 2*(3*x.^2-x+1);
P = @(x) x+x.^3/3; %antiderivative of p(x)=1+x^2

fprintf(['  h           error                error/h^2               h_max    ' ...
    '       error               error/h^2\n'])
for m = [9 99 999 9999]
    x = linspace(0,1,m+2)';
    h = x(2)-x(1);
    x2 = (((0:m+1)./(m+1)).^2)'; %nonuniform grid
    hmax = max(x2(2:m+2)-x2(1:m+1));
    x = x(2:m+1); x2 = x2(2:m+1);
    [u,errvec] = solve_fem(x,P,f,ufunc);
    [u2,errvec2] = solve_fem(x2,P,f,ufunc);
    err = norm(errvec,"inf"); err2 = norm(errvec2,"inf");
    fprintf('%.0e %.16d %.16d     %.4e %.16d %.16d\n',h,err,err/h^2,hmax,err2,err2/h^2)
end

function [u,errvec] = solve_fem(x,P,func,utrue)

m=length(x);
h = [x;1]-[0;x]; %ith entry of h denotes x_i-x_{i-1}

A=zeros(m);
A(1,1:2) = [(P(x(1))-P(0))/h(1)^2+(P(x(2))-P(x(1)))/h(2)^2, -(P(x(2))-P(x(1)))/h(2)^2];
for i = 2:m-1
    A(i,i-1:i+1) = [-(P(x(i))-P(x(i-1)))/h(i)^2, (P(x(i))-P(x(i-1)))/h(i)^2+...
        (P(x(i+1))-P(x(i)))/h(i+1)^2, -(P(x(i+1))-P(x(i)))/h(i+1)^2];
end
A(m,m-1:m) = [-(P(x(m))-P(x(m-1)))/h(m)^2, (P(x(m))-P(x(m-1)))/h(m)^2+...
    (P(1)-P(x(m)))/h(m+1)^2];

xmid = [x(1)/2; (x(1:m-1)+x(2:m))/2; (x(m)+1)/2]; %midpoints of subintervals
f = zeros(m,1);
f(1) = func(xmid(1))*(xmid(1)-0)+func(xmid(2))*(x(2)-xmid(2));
for i = 2:m-1
    f(i) = func(xmid(i))*(xmid(i)-x(i-1))+func(xmid(i+1))*(x(i+1)-xmid(i+1));
end
f(m) = func(xmid(m))*(xmid(m)-x(m-1))+func(xmid(m+1))*(1-xmid(m+1));

u=A\f;
errvec = u-utrue(x);
end