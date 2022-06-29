%This script solves problem 15.4 from ATAP
clear
i = 1; %vector count
for n = 3:3:100
    s = chebpts(n+2);
    s = s(2:n+1); %remove endpoints
    [~,Lconst] = lebesgue(s);
    Lvec(i) = Lconst;
    i = i+1;
end

plot(3:3:100,Lvec,'-ro')
xlabel('n')
ylabel('Lebesgue constant')
title('Lebesgue constant for Chebyshev points without endpoints')
saveas(gcf,'15-4','epsc')