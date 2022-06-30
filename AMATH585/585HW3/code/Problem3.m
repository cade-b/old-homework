clear
utrue=@(x) (1-x).^2;
ftrue=@(x) 2*(3*x.^2-2*x+1);
c = @(x) 1+x.^2;

m=9;
[u_coarse,h,errvec]=rodFD(c,utrue,ftrue,m);
errinf=norm(errvec,'inf');
fprintf('%.e       %.16e\n',h,errinf)

m=19;
[u_fine,h,errvec]=rodFD(c,utrue,ftrue,m);
errinf=norm(errvec,'inf');
fprintf('%.e       %.16e\n',h,errinf)

u_rich=(4*u_fine(1:2:end)-u_coarse)/3;
u_true=utrue(linspace(0,1,length(u_coarse)))';
errinf=norm(u_true-u_rich,'inf');
fprintf('Richardson       %.16e\n',errinf)


function [u,h,errvec] = rodFD(c,utrue,ftrue,m)
h=1/(m+1);
x = linspace(0,1,m+2)';
x = x(2:end); %remove x_0=0 fo indexing

A=zeros(m+1);
A(1,1:2)=[-(c(x(1)-h/2)+c(x(1)+h/2)), c(x(1)+h/2)]/h^2; %equation i=1
A(m+1,m:m+1)=[c(x(m+1)-h/2)+c(x(m+1)+h/2), -(c(x(m+1)-h/2)+c(x(m+1)+h/2))]/h^2; %i=m+1
for i = 2:m %interior
    A(i,i-1:i+1)=[c(x(i)-h/2), -(c(x(i)-h/2)+c(x(i)+h/2)), c(x(i)+h/2)]/h^2;
end

f = ftrue(x);
f(1)=f(1)-c(x(1)-h/2)/h^2; %Dirichlet BC

u=A\f;
u=[1;u]; %append u_0=1 to solution vector
x = [0;x]; %reappend x_0=0
uexact=utrue(x);
errvec=u-uexact;
end