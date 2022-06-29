%This is a modification of p5.m that uses only one FFT and IFFT to produce
%the same plots

% Differentiation of a hat function and exp(sin(x)):
  N = 24; h = 2*pi/N; x = h*(1:N)';
  t = max(0,1-abs(x-pi)/2); u = exp(sin(x)); %t: hat function, u:exp(sin(x))
  uprime = cos(x).*u;
  v = t+1i*u; %combined function v
  v_hat = fft(v);
  w_hat = 1i*[0:N/2-1 0 -N/2+1:-1]' .* v_hat;
  w = ifft(w_hat); 
  w_t= real(w); w_u= imag(w); clf
  subplot(3,2,1), plot(x,t,'.-','markersize',13)
  axis([0 2*pi -.5 1.5]), grid on, title('function')
  subplot(3,2,2), plot(x,w_t,'.-','markersize',13)
  axis([0 2*pi -1 1]), grid on, title('spectral derivative')
  subplot(3,2,3), plot(x,u,'.-','markersize',13)
  axis([0 2*pi 0 3]), grid on
  subplot(3,2,4), plot(x,w_u,'.-','markersize',13)
  axis([0 2*pi -2 2]), grid on
  error = norm(w_u-uprime,inf);
  text(2.2,1.4,['max error = ' num2str(error)])
  
  saveas(gcf,'prob3-6','epsc')
