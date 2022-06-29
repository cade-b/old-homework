%
%  Solves the steady-state heat equation in a square with conductivity
%  c(x,y) = 1 + x^2 + y^2:
%
%     -d/dx( (1+x^2+y^2) du/dx ) - d/dy( (1+x^2+y^2) du/dy ) = f(x),   
%                                                       0 < x,y < 1
%     u(x,0) = u(x,1) = u(0,y) = u(1,y) = 0
%
%  Uses a centered finite difference method.

%  Set up grid.

n = input(' Enter number of subintervals in each direction: ');
h = 1/n;
N = (n-1)^2;

%  Form block tridiagonal finite difference matrix A and right-hand side 
%  vector b.

A=sparse(zeros(N,N));
b = ones(N,1);         % Use right-hand side vector of all 1's.

%  Loop over grid points in y direction.

for j=1:n-1,
  yj = j*h;
  yjph = yj+h/2;  yjmh = yj-h/2;

%    Loop over grid points in x direction.

  for i=1:n-1,
    xi = i*h;
    xiph = xi+h/2;  ximh = xi-h/2;
    aiphj = 1 + xiph^2 + yj^2;
    aimhj = 1 + ximh^2 + yj^2;
    aijph = 1 + xi^2 + yjph^2;
    aijmh = 1 + xi^2 + yjmh^2;
    k = (j-1)*(n-1) + i;
    A(k,k) = aiphj+aimhj+aijph+aijmh;
    if i > 1, A(k,k-1) = -aimhj; end;
    if i < n-1, A(k,k+1) = -aiphj; end;
    if j > 1, A(k,k-(n-1)) = -aijmh; end;
    if j < n-1, A(k,k+(n-1)) = -aijph; end;
  end;
end;
A = (1/h^2)*A;   % Remember to multiply A by (1/h^2).

% Solve linear system.

u_comp = A\b;

x0 = zeros(length(b),1); %initial guess
tol = 1e-8; maxiter = 1000;

D = diag(diag(A));
[x_J,iter,residvec_J] = simpleiter(A,b,D,x0,tol,maxiter);
figure(1)
semilogy(0:iter-1,residvec_J,'LineWidth',1.6)
hold on
fprintf('Jacobi takes %i iterations.\n',iter)

L = tril(A);
[x_G,iter,residvec_G] = simpleiter(A,b,L,x0,tol,maxiter);
semilogy(0:iter-1,residvec_G,'--','LineWidth',1.6)
fprintf('Gauss-Seidel takes %i iterations.\n',iter)

omega1=1.4;
M1 = D/omega1+tril(A,-1);
[x_S1,iter,residvec_S1] = simpleiter(A,b,M1,x0,tol,maxiter);
semilogy(0:length(residvec_S1)-1,residvec_S1,':','LineWidth',1.6)
fprintf('SOR with omega=%d takes %i iterations.\n',omega1, iter)

omega2=1.6;
M2 = D/omega2+tril(A,-1);
[x_S2,iter,residvec_S2] = simpleiter(A,b,M2,x0,tol,maxiter);
semilogy(0:length(residvec_S2)-1,residvec_S2,':','LineWidth',1.6)
fprintf('SOR with omega=%d takes %i iterations.\n',omega2, iter)

omega3=1.9;
M3 = D/omega3+tril(A,-1);
[x_S3,iter,residvec_S3] = simpleiter(A,b,M3,x0,tol,maxiter);
semilogy(0:length(residvec_S3)-1,residvec_S3,':','LineWidth',1.6)
fprintf('SOR with omega=%d takes %i iterations.\n',omega3, iter)

% nonpreconditioned cg
[x_cg,~,~,iter,residvec_cg] = pcg(A,b,tol,maxiter,eye(size(A)));
semilogy(0:length(residvec_cg)-1,residvec_cg,'-.','LineWidth',1.6)
fprintf('CG takes %i iterations.\n',iter)

L = ichol(A); %incomplete Cholesky
[x_icg,~,~,iter,residvec_icg] = pcg(A,b,tol,maxiter,L,L');
semilogy(0:length(residvec_icg)-1,residvec_icg,'-.','LineWidth',1.6)
fprintf('ICCG takes %i iterations.\n',iter)
xlabel('iteraton'); ylabel('$\|b-Au^{(k)}\|/\|b\|$','Interpreter','latex')
legend('Jacobi','Gauss-Seidel',strcat('SOR with \omega=',num2str(omega1)) ...
    ,strcat('SOR with \omega=',num2str(omega2)),strcat('SOR with \omega=',num2str(omega3)) ...
    ,'unpreconditioned CG','ICCG')
title(strcat('Convergence of different methods with h=',num2str(h)))
hold off

saveas(gcf,strcat('h=',num2str(h),'.eps'),'epsc')