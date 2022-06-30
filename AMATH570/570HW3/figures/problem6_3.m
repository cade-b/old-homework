%This script produces a convergence plot for the Chebyshev approximations
%of a certain "difficult function"
clear
for i=1:13
    n = 2^(i+1)-1;
    nvec(i)=n;
    pn = chebfun(@(x) sin(1/x)*sin(1/sin(1/x)),[0.07,0.4],n);
    x = linspace(0.07,0.4,2*n);
    f = sin(1./x).*sin(1./sin(1./x));
    err(i)= norm(f-pn(x),inf);
end
figure
loglog(nvec,err,'-x')
hold on
loglog(nvec,nvec.^(-1/3))
legend('max errors','O(n^{-1/3})')
xlabel('n')
ylabel('max error')
title('Error of Chebyshev interpolant to sin(1/x)*sin(1/sin(1/x))')
hold off
saveas(gcf,'prob6-3','epsc')