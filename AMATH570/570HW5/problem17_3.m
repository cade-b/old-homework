%This script produces plots for parts a and b of exercise 17.8 in ATAP
clear
x = chebfun('x');
for k=0:5
    P{k+1}=x^k;
end
A = [P{1} P{2} P{3} P{4} P{5} P{6}]; %quasimatrix
[Q,R] = qr(A);
figure(1)
plot(Q)
xlabel('x')
ylabel('f(x)')
title('Columns of unmodified Q')
saveas(gcf,'17_1','epsc')
figure(2)
L_a = legpoly(0:5,'norm');
plot(L_a)
xlabel('x')
ylabel('f(x)')
title('Legendre polynomials normalized by (17.3)')
saveas(gcf,'17_2','epsc')
figure(3)
plot(abs(Q-L_a)); %error plot
xlabel('x')
ylabel('|Q-P_j|')
title('Difference between Q and Legendre polynomials')
saveas(gcf,'17_3','epsc')

%part b
for k=1:6
    Qk = Q(:,k);
    Q(:,k)=Q(:,k)/Qk(1); %normalize Q
    R(k,:)=R(k,:)*Qk(1); %adjust rows of R
end
Anew = Q*R;
figure(4)
plot(abs(A-Anew)) %comparison to original A
xlabel('x')
ylabel('|A-Q*R|')
title('Difference between original A and modified A')
saveas(gcf,'17_4','epsc')
figure(5)
plot(Q)
xlabel('x')
ylabel('f(x)')
title('Columns of modified Q')
saveas(gcf,'17_5','epsc')
figure(6)
L_b=legpoly(0:5);
plot(L_b)
xlabel('x')
ylabel('f(x)')
title('Legendre polynomials normalized by P_i(1)=1')
saveas(gcf,'17_6','epsc')
figure(7)
plot(abs(Q-L_b)); %error plot
xlabel('x')
ylabel('|Q-P_j|')
title('Difference between Q and Legendre polynomials')
saveas(gcf,'17_7','epsc')

