%This script determines the coefficient of T_5 in the Chebyshev expansion of arctan(x) on [−1, 1]
clear
x = chebfun('x');
f = atan(x);
coeffs=chebcoeffs(f);
fprintf('The coefficient of T_5 in the Chebyshev expansion of arctan(x) on [−1, 1] is approximately %d\n',coeffs(6));