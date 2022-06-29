clear

func = @(x,y) x.^2+y.^2;
lpfunc = @(x,y) 4*ones(length(x),length(y));

pow = 10; %power of 2 we use as true solution
m = 2^pow-1; h=1/(m+1);

x = linspace(0,1,m+2);
y = linspace(0,1,m+2);
[X,Y] = meshgrid(x,y);
x = x(2:end-1); y = y(2:end-1); %eliminate boundary


a = 1; %boundary condition

utrue = nineptpois(func,lpfunc,a,x,y); 

fprintf(['   h          L2 error          L2 ratio      ' ...
    '   Inf error        Inf ratio\n'])
for i = 2:6
m = 2^i-1; h=1/(m+1);
x = linspace(0,1,m+2);
y = linspace(0,1,m+2);
[X,Y] = meshgrid(x,y);
x = x(2:end-1); y = y(2:end-1); %eliminate boundary
u = nineptpois(func,lpfunc,a,x,y); 

utrue_adj = utrue(2^(pow-i):2^(pow-i):end,2^(pow-i):2^(pow-i):end);
err = u-utrue_adj;
err_L2 = h*norm(err(:),2);
err_inf = norm(err(:),'inf');
constvec = err_inf/h^2;
constvec2 = err_L2/h^3;
fprintf('%1.4e %.10e %.10e %.10e %.10e\n',h,err_L2, constvec2, err_inf, constvec)
end


function [u,uBC] = nineptpois(func,lpfunc,a,x,y)
m = length(x); h = 1/(m+1);
%From pg. 68 of course text
I = eye(m);
e = ones(m,1);
T = spdiags([4*e -20*e 4*e],[-1 0 1],m,m);
S1 = spdiags([e e],[-1 1],m,m);
S2 = spdiags([e 4*e e],[-1 0 1],m,m);
A = (kron(I,T) + kron(S1,S2))/(6*h^2);

f = func(x',y); 
%implement boundary condition
f(1,:) = f(1,:) - a/h^2;
f(m,:) = f(m,:) - a/h^2;
f(:,1) = f(:,1) - a/h^2;
f(:,m) = f(:,m) - a/h^2;
f(1,1) = f(1,1) + a/(6*h^2);
f(1,m) = f(1,m) + a/(6*h^2);
f(m,1) = f(m,1) + a/(6*h^2);
f(m,m) = f(m,m) + a/(6*h^2);

lpf = lpfunc(x',y);
f = f + lpf*h^2/12;

f = f(:);

u = A\f;
u = reshape(u,m,m);
uBC = a*ones(m+2);
uBC(2:end-1,2:end-1)=u; %add back boundary condition

end