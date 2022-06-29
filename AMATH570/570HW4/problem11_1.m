%This code verifies the Hermite integral formula in a particular instance
%by computing contour integrals in chebfun (problem 11.1 in ATAP)
clear
s = chebfun(@(s) s,[0 2*pi]); 
l = chebfun(@(x) (x-0.5)*(x-1)*(x+1)); %node polynomial
z_1 = chebfun('1.5*exp(1i*s)',[0 2*pi]); %|x|=3/2
z_2 = chebfun('3*exp(1i*s)',[0 2*pi]); %|x|=3
f = chebfun(@(x) (x+1)*(x-0.5)*(x-1)*exp(x)+11/6+x/2-x^2/3); 

x=2; %evaluate at x=2
%|x|=3/2 contour
inner_term_1 = f(z_1)*(l(z_1)-l(x))/(l(z_1)*(z_1-x)); %term inside integral
output1 = sum(inner_term_1*diff(z_1))/(2i*pi);
fprintf('For the contour |x|=3/2, p(%i)=%d\n',x,output1)

%|x|=3 contour
inner_term_2 = f(z_2)*(l(z_2)-l(x))/(l(z_2)*(z_2-x));
output2 = sum(inner_term_2*diff(z_2))/(2i*pi);
fprintf('For the contour |x|=3, p(%i)=%d\n',x,output2)

%part c
p = @(x) sum(f(z_1)*(l(z_1)-l(x))/(l(z_1)*(z_1-x))*diff(z_1))/(2i*pi);
pcheb = chebfun(@(x) p(x));
coeff = poly(pcheb);
disp('The coefficients that poly gives are:')
disp(coeff)


