%This script plots various upper bounds for the convergence of an entire function 
clear
f = @(x) exp(-x^2);
fexact = chebfun(f); 
for n = 0:length(fexact)-1
    pn = chebfun(f,n);
    err(n+1) = norm(pn-fexact,inf);
end
figure
semilogy(0:length(fexact)-1,err,'-x')
xlabel('n')
ylabel('max error')
title('Convergence of exp(−x^2) with upper bounds')
hold on

for rho = [1.1, 1.2, 1.4, 2, 3, 5, 8]
    M = exp(((rho-1/rho)/2)^2);
    upbound = 4*M*rho.^(-(0:length(fexact)-1))/(rho-1);
    semilogy(0:length(fexact)-1,upbound)
end
legend('max error','\rho=1.1','\rho=1.2','\rho=1.4','\rho=2','\rho=3','\rho=5','\rho=8')
hold off
saveas(gcf,'prob8-3','epsc')

