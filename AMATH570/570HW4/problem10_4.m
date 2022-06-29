%This script solves problem 10.4 in SMM
clear
N = 16; dt = 0.00001; %chosen for stability
tmax = 5; %arbitrarily chosen
[D,x] = cheb(N);
D2 = D^2; D2 = D2(2:N,2:N); %get operator
nmax = round(tmax/dt);
udata = zeros(N+1,nmax);
tvec = dt*(0:(nmax-1)); %keep track of values of t
for n=2:nmax
    %RK4 scheme
    k1 = dt*(D2*udata(2:N,n-1)+exp(udata(2:N,n-1)));
    k2 = dt*(D2*(udata(2:N,n-1)+k1/2)+exp(udata(2:N,n-1)+k1/2));
    k3 = dt*((D2*(udata(2:N,n-1)+k2/2)+exp(udata(2:N,n-1)+k2/2)));
    k4 = dt*((D2*(udata(2:N,n-1)+k3)+exp(udata(2:N,n-1)+k3)));
    udata(2:N,n)=udata(2:N,n-1)+(k1 + 2*k2 + 2*k3 + k4)/6;
end

t35 = 3.5/dt+1;
fprintf('u(0,3.5) is approximately %d\n',udata(N/2+1,t35))

k = find(udata(N/2+1,:)>=5,1);
fprintf('u(0,t)=5 at time t=%d\n',tvec(k))

