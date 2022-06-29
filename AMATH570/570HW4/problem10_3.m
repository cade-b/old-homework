%This script solves problem 10.3 in SMM
clear
for nu = 7:8
N=50;
[D,x] = cheb(N);
u = exp(-60*(x - 1/2).^2);
tmax = 1; 
dt=nu*N^(-2); 
A = -D(1:N,1:N); %differential operator for this problem
nmax = round(tmax/dt);
udata = zeros(N+1,nmax); %store solution
udata(:,1) = u;
udata(:,2)=exp(-60*(x + dt - 1/2).^2); %using exact solution
udata(:,3)=exp(-60*(x + 2*dt - 1/2).^2); 
udata(end,:)=0; %boundary condition
for n = 4:nmax
    %AB3 scheme
    udata(1:N,n)=udata(1:N,n-1)+(dt/12)*(23*A*udata(1:N,n-1)-16*A*udata(1:N,n-2)+5*A*udata(1:N,n-3));
    udata(N+1,n) = 0; %boundary condition
end

figure
plot(x,udata(:,end))
xlabel('x')
ylabel('u(x,1)')
title('Computed solution at t=1')
saveas(gcf,sprintf('sol%d',nu),'epsc')

%second plot
figure
z = exp(1i*pi*(0:200)/100); r = z-1;
s = (23-16./z+5./z.^2)/12; plot(r./s) %stability region
hold on
plot(eig(A)*dt,'.') %lambda*dt
x = -1:.02:0.25; y = -0.8:.02:0.8; [xx,yy] = meshgrid(x,y); zz = xx+1i*yy;
I = eye(N); sigmin = zeros(length(y),length(x));
for j = 1:length(x)
    for i = 1:length(y), sigmin(i,j) = min(svd(zz(i,j)*I-A*dt)); end
end
contour(x,y,sigmin,10.^(-6:1:-2)), colormap(1e-6*[1 1 1]);
grid on
hold off
title('Stability region, \lambda*dt, and \epsilon-pseudospectra')
saveas(gcf,sprintf('reg%d',nu),'epsc')
end
