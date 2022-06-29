%This script solves problem 2 on hw 
clear
x = chebfun('x');
f_a = 1/(1+25*x^2);
f_b = 1/(1+500*x^2);

i = 1;
nvec = 2:50;
for n = nvec
    %Chebyshev points
    scheb = chebpts(n+2);
    [~,Lconst_cheb] = lebesgue(scheb);
    [~,err_a] = minimax(f_a,n+2);
    bound_cheb_a(i) = (Lconst_cheb+1)*err_a;
    
    [~,err_b] = minimax(f_b,n+2);
    bound_cheb_b(i) = (Lconst_cheb+1)*err_b;
    
    p_a_cheb = chebfun(f_a,n+3);
    err_a_cheb(i) = norm(f_a-p_a_cheb,inf);
    
    p_b_cheb = chebfun(f_b,n+3);
    err_b_cheb(i) = norm(f_b-p_b_cheb,inf);
    
    %Equispaced points
    sequi = linspace(-1,1,n+2)';
    [~,Lconst_equi] = lebesgue(sequi);
    bound_equi_a(i) = (Lconst_equi+1)*err_a;
    
    bound_equi_b(i) = (Lconst_equi+1)*err_b; 
    
    p_a_equi = interp1(sequi,f_a(sequi),domain(-1,1));
    err_a_equi(i) = norm(f_a-p_a_equi,inf);
    
    p_b_equi = interp1(sequi,f_b(sequi),domain(-1,1));
    err_b_equi(i) = norm(f_b-p_b_equi,inf);
    
    i = i+1;
end

figure(1)
semilogy(nvec,bound_cheb_a)
hold on
semilogy(nvec,err_a_cheb)

semilogy(nvec,bound_equi_a)
semilogy(nvec,err_a_equi)
xlabel('n')
ylabel('error')
legend('bound for Chebyshev points','Chebyshev points','bound for equispaced points','equispaced points')
title('Bound 15.5 for f(x)=1/(1+25*x^2)')
hold off
saveas(gcf,'2a','epsc')

figure(2)
semilogy(nvec,bound_cheb_b)
hold on
semilogy(nvec,err_b_cheb)

semilogy(nvec,bound_equi_b)
semilogy(nvec,err_b_equi)
xlabel('n')
ylabel('error')
legend('bound for Chebyshev points','Chebyshev points','bound for equispaced points','equispaced points')
title('Bound 15.5 for f(x)=1/(1+500*x^2)')
hold off
saveas(gcf,'2b','epsc')