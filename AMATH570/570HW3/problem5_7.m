%this script chooses interpolation points via a greedy algorithm 
clear
f = @(x) abs(x);
greed(f)

function greed(func)
    x = chebfun('x');
    f = func(x);
    d = domain(-1,1); 
    p = 0;
    s = [];
    for k = 0:25
        [~,maxval]=max(abs(f-p));
        s = [s, maxval];
        y = f(s);
        p = interp1(s,y,d);
        err = f-p;
        plot(err)
        plot(s,err(s),'.k')
        hold on
    end
    fprintf('Final gridpoints are:')
    disp(sort(s))
    xlabel('x')
    ylabel('(f-p)(x)')
    title('Error plots for greedy algorithm')
    axis([-1,1,-inf,inf])
    hold off
    saveas(gcf,'prob5-7a','epsc')
    
    figure(2)
    plot(err)
    hold on
    plot(s,err(s),'.k')
    xlabel('x')
    ylabel('(f-p)(x)')
    title('Final error plot for greedy algorithm')
    hold off
    saveas(gcf,'prob5-7b','epsc')
end

