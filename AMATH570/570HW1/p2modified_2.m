% This is a modified version of p2.m to take u(x) = exp(sin(x)*abs(sin(x)))

% For various N (even), set up grid as before:
  clf, subplot('position',[.1 .4 .8 .5])
  for N = 2:2:100;
    h = 2*pi/N;
    x = -pi + (1:N)'*h;
    u = exp(sin(x).*abs(sin(x))); uprime = 2*cos(x).*abs(sin(x)).*u;

    % Construct spectral differentiation matrix:
    column = [0 .5*(-1).^(1:N-1).*cot((1:N-1)*h/2)];
    D = toeplitz(column,column([1 N:-1:2]));

    % Plot max(abs(D*u-uprime)):
    error = norm(D*u-uprime,inf);
    loglog(N,error,'.','markersize',15), hold on
  end
  grid on, xlabel N, ylabel error
  semilogy(2:2:100,(2:2:100).^(-1),'--') 
  text(50,5e-4,'N^{-1}','fontsize',18)
  title('Convergence of spectral differentiation')
  saveas(gcf,'p2b','epsc')
