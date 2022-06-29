clear

T=2*pi;
m=2001;
t=linspace(0,T,m+2)';
h = t(2)-t(1);
alpha=0.7; beta=0.7;

theta=0.7+t.*(t-2*pi).^2;


theta = theta(2:end-1); %remove boundary
tol=1e-12; itermax=1000;
[theta,~,~]=NewtonSolve(theta,alpha,beta,h,tol,itermax);
figure(1)
plot(t,theta)
plot(t,theta)
xlabel('t'); ylabel('\theta(t)')
title('New solution for given BCs')
saveas(gcf,'hw3p1a','epsc')

i=1;
fprintf('T           max(theta)\n')

theta=0.7+sin(t/2);
theta = theta(2:end-1); %remove boundary
[theta,~,~]=NewtonSolve(theta,alpha,beta,h,tol,itermax);
fprintf('%i      %16.16d\n',T,max(theta))
figure(2)
for T=[20,50,100,200]
    t=linspace(0,T,m+2)';
    h = t(2)-t(1);
    theta = theta(2:end-1); %remove boundary, reuse old theta as new guess
    [theta,~,~]=NewtonSolve(theta,alpha,beta,h,tol,itermax);
    subplot(2,2,i)
    plot(t,theta)
    title(['T=',num2str(T)])
    xlabel('t')
    ylabel('\theta(t)')
    fprintf('%i      %16.16d\n',T,max(theta))
    i=i+1;   
end
saveas(gcf,'hw3p1b','epsc')

function G=buildG(theta,alpha,beta,h)

G=zeros(length(theta),1);
theta = [alpha; theta; beta]; %include BCs for computation 
for i=1:length(theta)-2
    G(i)=(theta(i)-2*theta(i+1)+theta(i+2))/h^2+sin(theta(i+1));
end

end


function J=buildJ(theta,h)

maindiag=-2+cos(theta)*h^2;
J=spdiags([ones(length(theta),1),maindiag,ones(length(theta),1)],-1:1,...
    length(theta),length(theta));
J=J/h^2;

end


function [theta,iter,dnormvec]=NewtonSolve(initial_guess,alpha,beta,h,tol,itermax)

theta = initial_guess;
deltanorm=Inf;
iter=0;
while deltanorm>tol && iter<itermax
    G=buildG(theta,alpha,beta,h);
    J=buildJ(theta,h);
    delta=-J\G;
    theta=theta+delta;
    deltanorm=norm(delta,'inf');
    iter=iter+1;
    dnormvec(iter)=deltanorm;
end

theta=[alpha;theta;beta]; %add BCs

end