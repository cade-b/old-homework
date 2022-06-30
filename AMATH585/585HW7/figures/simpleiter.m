function [x,iter,residvec,X] = simpleiter(A,b,M,x0,tol,maxiter)
x = x0; r = b-A*x0;
z = M\r;
errnorm = norm(r);
iter = 0;
nb = norm(b); %save norm(b)
while errnorm/nb>tol && iter<maxiter
    x = x + z;
    r = b - A*x;
    z = M\r;
    iter = iter + 1;
    errnorm = norm(r);
    residvec(iter) = errnorm/nb; %relative residual
    X(:,iter) = x;
end
end