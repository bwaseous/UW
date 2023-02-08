function plot_rhr(f, a, b, dx)
    xplot = a:(dx/100):b;
    x = a:dx:b;
    n = length(x);
    y = f(x);
    
    plot(xplot, f(xplot), 'k')
    hold on
    
    for k = 1:n-1
        x_rect = [x(k) x(k) x(k+1) x(k+1) x(k)];
        y_rect = [0 f(x(k+1)) f(x(k+1)) 0 0];
        plot(x_rect, y_rect, 'b')
    end
    hold off
end