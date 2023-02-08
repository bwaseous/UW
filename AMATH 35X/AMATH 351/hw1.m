[x, y] = meshgrid(-10:0.4:10, -20:0.4 :20); 
s = y - x.^2;    % for slope function f = x^2*cos(y)
L = sqrt(1 + s.^2);

close all

quiver(x,y, 1./L, s./L, 0.5), axis tight
xlabel('x');
ylabel('y');
title('Direction field for y'' = (x^2 - 1)(y^2 - 1)');

hold on
xspan=-10:0.5:10;
plot(xspan,-1*ones(size(xspan)),'b-')
plot(xspan,1*ones(size(xspan)),'b-')
%plot(xspan,-pi*ones(size(xspan)),'b-')
%plot(xspan,-2*pi*ones(size(xspan)),'b-')

tspan = [0,3];

[t,y1] = ode45(@(t,y1) (y1 - t.^2), tspan, 3);
plot(t,y1,'k-')

[t,y2] = ode45(@(t,y2) (y2 - t.^2), tspan, -1);
plot(t,y2,'k-')

%[t,y3] = ode45(@(t,y3) t^2*sin(y3), tspan, -3);
%plot(t,y3,'k-')

ax = gca;
ax.FontSize = 12;
ax.FontName = 'Times';

hold off