%This script finds the total variation of sin(100*x)/(1+x^2)
clear
tv = @(f) norm(diff(f),1); %tv function from text
x = chebfun('x');
f= sin(100*x)/(1+x^2);
fprintf('The total variation of f(x) is %d.\n',tv(f))
