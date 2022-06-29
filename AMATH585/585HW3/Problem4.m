clear 

figure(1)
for m=[19 39 79 159] 
    a=0; b=1; alpha=-1; beta=1.5; epsilon=0.01;
    x=linspace(a,b,m+2)';
    h = x(2)-x(1);
    tol=1e-12; itermax=100;

    initial_guess=buildinitialguess(x(2:end-1),a,b,alpha,beta,epsilon);
    [u,iter,dnormvec]=NewtonSolve(initial_guess,alpha,beta,epsilon,h,tol,itermax);
    plot(x,u)
    hold on
end
xlabel('x')
ylabel('u(x)')
legend('h=1/20','h=1/40','h=1/80','h=1/160')
title('Solution of difference equations for different h')
hold off
saveas(gcf,'hw3p4','epsc')

function u=buildinitialguess(x,a,b,alpha,beta,epsilon) %function from 2.105

w0=(a-b+beta-alpha)/2;
xbar=(a+b-alpha-beta)/2;
u=x-xbar+w0*tanh(w0*(x-xbar)/(2*epsilon));

end



function G=buildG(u,alpha,beta,h,epsilon)

G=zeros(length(u),1);
u = [alpha; u; beta];
for i=1:length(u)-2
    G(i)=epsilon*(u(i)-2*u(i+1)+u(i+2))/h^2+u(i+1)*((u(i+2)-u(i))/(2*h)-1);
end

end



function J=buildJ(u,alpha,beta,h,epsilon)

u = [alpha; u; beta];

for i=1:length(u)-2
    Jdiag(i)=-2*epsilon/h^2+(u(i+2)-u(i))/(2*h)-1;
end
Jsuperdiag=epsilon/h^2+u(1:end-2)/(2*h);
Jsubdiag=epsilon/h^2-u(3:end)/(2*h);

J=spdiags([Jsubdiag,Jdiag',Jsuperdiag],-1:1,length(u)-2,length(u)-2);

end


function [u,iter,dnormvec]=NewtonSolve(initial_guess,alpha,beta,epsilon,h,tol,itermax)

u = initial_guess;
deltanorm=Inf;
iter=0;
while deltanorm>tol && iter<itermax
    G=buildG(u,alpha,beta,h,epsilon);
    J=buildJ(u,alpha,beta,h,epsilon);
    delta=-J\G;
    u=u+delta;
    deltanorm=norm(delta,'inf');
    iter=iter+1;
    dnormvec(iter)=deltanorm;
end

u=[alpha;u;beta]; %add BCs

end