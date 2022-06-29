%This script plots the max error for the Chebyshev approximations of exp(x)
%and Runge's function for various vaules of n
clear
x = chebfun('x');
f_a = exp(x);
f_b = 1/(1+25*x^2);
for n = 1:200
    p_a = chebfun(f_a,n+1);
    p_b = chebfun(f_b,n+1);
    error_a(n) = norm(f_a-p_a,inf);
    error_b(n) = norm(f_b-p_b,inf);
end
figure
semilogy(1:200,error_a)
hold on
semilogy(1:200,error_b)
xlabel('n')
ylabel('||f−p_n||')	
legend('e^x','1/(1+25x^2)')
title('Max error of Chebyshev interpolants')
saveas(gcf,'prob2-5','epsc')