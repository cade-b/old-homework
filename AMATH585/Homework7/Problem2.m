clear

phi = @(x) 20*pi*x.^3;
dphi = @(x) 60*pi*x.^2;
ddphi = @(x) 120*pi*x;
a = 1/2;
func = @(x) -20+a*ddphi(x).*cos(phi(x))-a*(dphi(x).^2).*sin(phi(x));
ufunc = @(x) 1+12*x-10*x.^2+a*sin(phi(x));

m = 255; h = 1/(m+1); 
x = linspace(0,1,m+2)';
x = x(2:end-1);
e = ones(m,1);
A=-spdiags([e -2*e e],[-1,0,1],m,m)/h^2;
f = func(x);
f(1) = f(1)-1/h^2;
f(end) = f(end)-3/h^2;
f = -f; %need to negate both sides to ensure positive definiteness for cg

utrue = A\f; %true solution to linear system

u0 = 1+2*x;
tol = -1; %ensure tolerance isn't hit so we run all 20 iterations
maxiter = 20;

M = tril(A); %Gauss-Seidel
[~,~,residvec,U_gs] = simpleiter(A,f,M,u0,tol,maxiter);

for iter = 1:maxiter
    U_cg(:,iter) = pcg(A,f,tol,iter,eye(m),eye(m),u0);
end

fprintf(['iterations   G-S L2 error   G-S infinity error   CG L2 error' ...
    '  CG Infinity Error\n'])
err = u0-utrue;
err2 = sqrt(h)*norm(err); 
errinf = norm(err,'inf');
fprintf('    0      %.10e  %.10e %.10e %.10e\n',err2,errinf,err2,errinf)
figure()
plot(x,u0,x,utrue)
xlabel('x'); ylabel('u(x)');
title('Solution at iteration 0 for both methods')
legend('approximate solution','true solution')
saveas(gcf,'initial.eps','epsc')

figure()
plot(x,err)
xlabel('x'); ylabel('u(x)');
title('Error at iteration 0 for both methods')
saveas(gcf,'initialerr.eps','epsc')

for i = [5 10 20]
    u_gs = U_gs(:,i);
    u_cg = U_cg(:,i);
    err_gs = u_gs-utrue;
    err2_gs = sqrt(h)*norm(err_gs);
    errinf_gs = norm(err_gs,'inf');
    err_cg = u_cg-utrue;
    err2_cg = sqrt(h)*norm(err_cg);
    errinf_cg = norm(err_cg,'inf');
    fprintf('    %i      %.10e  %.10e %.10e %.10e\n',i, ...
        err2_gs,errinf_gs,err2_cg,errinf_cg)

    figure()
    plot(x,u_gs,x,utrue)
    xlabel('x'); ylabel('u(x)');
    title(['Solution at iteration ',num2str(i),' for G-S'])
    legend('approximate solution','true solution')
    saveas(gcf,strcat('GS_i=',num2str(i),'.eps'),'epsc')

    figure()
    plot(x,err_gs)
    xlabel('x'); ylabel('u(x)');
    title(['Error at iteration ',num2str(i),' for G-S'])
    saveas(gcf,strcat('GSerr_i=',num2str(i),'.eps'),'epsc')

    figure()
    plot(x,u_cg,x,utrue)
    xlabel('x'); ylabel('Error');
    title(['Solution at iteration ',num2str(i),' for CG'])
    legend('approximate solution','true solution')
    saveas(gcf,strcat('CG_i=',num2str(i),'.eps'),'epsc')

    figure()
    plot(x,err_cg)
    xlabel('x'); ylabel('Error');
    title(['Error at iteration ',num2str(i),' for CG'])
    saveas(gcf,strcat('CGerr_i=',num2str(i),'.eps'),'epsc')
end
