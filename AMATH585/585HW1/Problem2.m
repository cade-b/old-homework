clear;

u=@(x) sin(x);
x = pi/6;
upptrue = -sin(x); %true u''(x) value
fprintf(' h       phi_0(h)       Error\n')

h = [0.2 0.1 0.05];
upp = (u(x+h)+u(x-h)-2*u(x))./h.^2; %phi_0
err = upp-upptrue;

for i=1:length(h)
    fprintf('%0.2f   %d %d\n',h(i), upp(i), err(i))
end

%phi_1
fprintf(' h       phi_1(h)        Error\n')
R1 = (4*upp(2)-upp(1))/3; %combine 0.2 and 0.1
errR1=R1-upptrue;
R2 = (4*upp(3)-upp(2))/3; %combine 0.1 and 0.05
errR2=R2-upptrue;
fprintf('%0.2f  %d  %d\n',h(1), R1, errR1)
fprintf('%0.2f  %d  %d\n',h(2), R2, errR2)

R3 = (16*R2-R1)/15; %phi_2
errR3=R3-upptrue;
fprintf(' h       phi_2(h)        Error\n')
fprintf('%0.2f  %d  %d\n',h(1), R3, errR3)

