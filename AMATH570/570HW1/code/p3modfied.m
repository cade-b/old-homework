%This is a modified version of p3.m that produces an error plot
clear
xmax = 10; clf
h=2^(-6);
xx = -xmax-h/20:h/10:xmax+h/20; % plotting grid; set to account for smallest value of h for consistency
for j=1:4
    h = 2^(-j-2);
    hvec(j)=h;
    x = -xmax:h:xmax;          % computational grid
    v1 = (abs(x)<=3);          % square wave
    vv1 = (abs(xx)<=3);        % square wave on plotting grid
    v2 = max(0,1-abs(x)/3);    % hat function
    vv2 = max(0,1-abs(xx)/3);  % hat function on plotting grid
    p1 = zeros(size(xx));
    p2 = zeros(size(xx));
    for i = 1:length(x),
      p1 = p1 + v1(i)*sin(pi*(xx-x(i))/h)./(pi*(xx-x(i))/h);
      p2 = p2 + v2(i)*sin(pi*(xx-x(i))/h)./(pi*(xx-x(i))/h);
    end
    errvec1(j)=norm(vv1-p1,inf);
    errvec2(j)=norm(vv2-p2,inf);
end
figure(1)
loglog(hvec,errvec1,'rx-')
hold on
loglog(hvec,errvec2,'bo-')
loglog(hvec,hvec/2)
xlabel('h')
ylabel('maximum error')
legend('square wave','hat function','O(h)')
title('Maximum error in sinc function interpolants')
saveas(gcf,'p3','epsc')
