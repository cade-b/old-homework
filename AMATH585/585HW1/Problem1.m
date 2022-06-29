clear;

u=@(x) sin(x);
x = pi/6;
upptrue = -sin(x); %true u''(x) value
fprintf('  h       FD Quotient      Error\n')
for k=1:16
    h = 10^-k;
    upp = (u(x+h)+u(x-h)-2*u(x))/h^2;
    err = upp-upptrue; %error term
    fprintf('%.e   %e  %e\n',h, upp, err)
end



