%This script compares the results of different polynomial fits
clear
fprintf('            Chebfun           polyfit         polyfitA     condition number\n')
for i = 1:10
    %part a
    k=i*10;
    f = @(x) cos(k*x);
    xx = (-1:0.001:1)'; %fine grid
    x = chebfun('x');
    p = f(x);
    n = length(p);
    pval = polyval(polyfit(xx,f(xx),n),0);
    [d, H] = polyfitA(xx,f(xx),n);
    pvalA = polyvalA(d,H,0);
    
    %part b
    %Vandermonde matrix creation from polyfit function
    V = zeros(length(xx),n+1);
    V(:,n+1) = ones(length(xx),1);
    for j = n:-1:1
        V(:,j) = xx.*V(:,j+1);
    end
    condnum = cond(V);
    
    %print results in one table
    fprintf('k=%i: %16d %16d %16d %16e\n',k,p(0),pval,pvalA,condnum)
end