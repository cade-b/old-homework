clear

phi = @(x) 20*pi*x.^3;
dphi = @(x) 60*pi*x.^2;
ddphi = @(x) 120*pi*x;
a = 1/2;
func = @(x) -20+a*ddphi(x).*cos(phi(x))-a*(dphi(x).^2).*sin(phi(x));
ufunc = @(x) 1+12*x-10*x.^2+a*sin(phi(x))-(1+2*x);

tol = 1e-13; maxiter=1000;

for m = [19 49 99 999]
    h = 1/(m+1); %need m odd here
    m_c = (m-1)/2; h_c = 2*h;
    x = linspace(0,1,m+2)';
    x = x(2:end-1);
    x_c = x(2:2:end);
    f = func(x); %randomly generate RHS vector


    % build FD matrix
    e = ones(m,1);
    A=spdiags([e -2*e e],[-1,0,1],m,m)/h^2;
    e_c = ones(m_c,1);
    A_c=spdiags([e_c -2*e_c e_c],[-1,0,1],m_c,m_c)/h_c^2;
    utrue = A\f;

    omega = 2/3; %suggested weight from class
    M = diag(diag(A))/omega; %weighted Jacobi

    u0 = zeros(size(utrue));

    [u,errvec,iter] = vcycle(A,f,A_c,M,u0,tol,maxiter);
    fprintf('Weighted Jacobi converged in %i iterations for h=%.2d.\n',iter,h)

    M=tril(A);
    [u,errvec,iter] = vcycle(A,f,A_c,M,u0,tol,maxiter);
    fprintf('Gauss-Seidel converged in %i iterations for h=%.2d.\n',iter,h)
end

function [u,resvec,iter] = vcycle(A,f,A_c,M,u0,tol,maxiter)
m = length(A);

e = ones(m,1);
Icf=spdiags([e/2 e e/2],[-1,0,1],m,m);
Icf = Icf(:,2:2:end);
Ifc = Icf'/2; 

r = f-A*u0;
u = u0; iter=0;
nf = norm(f); %save norm f
resnorm=norm(r)/nf; 
while resnorm>tol && iter<maxiter
    %beginning smoothing step
    u = u+M\r;
    r = f-A*u;

    r_c = Ifc*r;
    z_c = A_c\r_c;
    z = Icf*z_c;
    u = u+z;
    r = f-A*u;

    %ending smoothing step
    u = u+M\r;
    r = f-A*u;
    resnorm = norm(r)/nf;
    iter = iter+1;
    resvec(iter) = resnorm;
end

end



