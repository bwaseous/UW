%%% Problem 1
%%% Implement the Bisection method as we did in the Week 5 Coding Lecture
t = (1:0.001:3);
c = 1.3*((exp(1)).^(t/-11) - (exp(1)).^(-4*t/3));
    %plot(t,c,'linewidth', 3);
    %a=1;b=3;
    %dc = 1.3*((-1/11)*(exp(1).^(t/-11)) + (4/3)*((exp(1)).^(-4*t/3)));

a = 1;
b = 3;

dc_a = 1.3*((-1/11)*(exp(1).^(a/-11)) + (4/3)*((exp(1)).^(-4*a/3)));
dc_b = 1.3*((-1/11)*(exp(1).^(b/-11)) + (4/3)*((exp(1)).^(-4*b/3)));
xmid = (a+b)/2;
dc_mid = 1.3*((-1/11)*(exp(1).^(xmid/-11)) + (4/3)*((exp(1)).^(-4*xmid/3)));

for i=1:length(t)
    if abs(dc_mid) < 1e-8
        break
    else
        if sign(dc_mid) == sign(dc_a)
            a = xmid;
        else
            b = xmid;
        end
    end
    xmid = (a+b)/2;
    dc_mid = 1.3*((-1/11)*(exp(1).^(xmid/-11)) + (4/3)*((exp(1)).^(-4*xmid/3)));
    dc_a = 1.3*((-1/11)*(exp(1).^(a/-11)) + (4/3)*((exp(1)).^(-4*a/3)));
    dc_b = 1.3*((-1/11)*(exp(1).^(b/-11)) + (4/3)*((exp(1)).^(-4*b/3)));
end

A1 = xmid;
A2 = 1.3*((exp(1)).^(xmid/-11) - (exp(1)).^(-4*xmid/3));
A3 = abs(dc_mid);

%%% Problem 2
%%% Implement Newton's method as we did in the Week 5 Coding Lecture

xold = 2;
for i=1:1e8
    dx = 2.*xold;
    ddx = 2;
    xnew = xold - (dx)./(ddx);
    T = abs(xnew - xold);
    xold = xnew;
    if T<1e-8
        break
    end
end

A4 = i;
A5 = xnew;


xold = 2;
for i=1:1e8
    dx = 500.*(xold.^499);
    ddx = (500*499).*(xold^498);
    xnew = xold - (dx)./(ddx);
    T = abs((500.*(xnew.^499)) - (500.*(xold.^499)));
    xold = xnew;
    if T<1e-8
        break
    end
end
    
A6 = i+1;
A7 = xnew;


xold = 2;
for i=1:1e8
    dx = 1000.*(xold.^999);
    ddx = (1000*999).*(xold^998);
    xnew = xold - (dx)./(ddx);
    T = abs((1000.*(xnew.^999)) - (1000.*(xold.^999)));
    xold = xnew;
    if T<1e-8
        break
    end
end
    
A8 = i+1;
A9 = xnew;