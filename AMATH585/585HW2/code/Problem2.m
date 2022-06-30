clear
h=1/3;
x = [0 1/3 2/3 1];
u = [1 77/81 76/81 1];

G0 = [1 2/3 1/3 0];
G1 = [0 1/3 2/3 1]; 
G1_3 = [0 -2/27 -1/27 0]/h;
G2_3 = [0 -1/27 -2/27 0]/h;

plot(x,u,'x-',x,G0,'x-',x,G1,'x-',x,G1_3,'x-',x,G2_3,'x-')
legend('solution approximation','G_0(x)','G_1(x)','G(x,1/3)','G(x,2/3)')
xlabel('x')
title('Approximate solution and associated Green''s functions')
saveas(gcf,'hw2p2','epsc')
