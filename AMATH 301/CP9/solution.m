%%% Problem 1
%%% Use ode45 to solve the Fitzhugh-Nagumo IVP
%%% For the maxima use the plot to narrow down the times you use to search
%%% for the maximum.

v0 = [1;0];
alpha = 0.7;
beta = 1;
tau = 12;
params = [alpha,beta,tau];
tspan = (0:0.5:100);

[T,V] = ode45(@(t,v) FHN(t,v,params),tspan,v0);

A1 = V(:,1);
%plot(tspan, A1)
[M, I] = max(A1(1:30,:));
A2 = T(I);
[M, I] = max(A1(80:100,:));
A3 = T(I+79);
A4 = 1/(A3-A2);

%%% Problem 2
%%% Use ode45 to solve the Chua equation
%%% You can tell something is chaotic if it is seemingly random
%%% If it looks like all solutions tend toward a point or a circle it is
%%% not chaotic.

x0 = [0.1;0.2;0.3];
alpha = 16;
beta = 30;
params = [alpha,beta];
tspan = (0:0.05:100);

[T,X] = ode45(@(t,x) Chuas(t,x,params),tspan,x0);

%x1 = X(:,1); y1 = X(:,2); z1 = X(:,3);
%plot3(x1,y1,z1)
A5 = 1;
A6 = X';

x0 = [0.1;0.2+10e-5;0.3];

[T,X] = ode45(@(t,x) Chuas(t,x,params),tspan,x0);

A7 = max(max(abs(X'-A6))); 

x0 = [0.1;0.2;0.3];
beta = 100;
params = [alpha,beta];

[T,X] = ode45(@(t,x) Chuas(t,x,params),tspan,x0);

%x1 = X(:,1); y1 = X(:,2); z1 = X(:,3);
%plot3(x1,y1,z1)

A8 = 0;

%%% Problem 3

%%% Part 1: Finite Differences
%%% Use finite differences to solve the BVP
%%% Be careful about the shape of the vectors, you may have to transpose to
%%% get the correct shape.  It's a good idea to print the solutions out to
%%% make sure the shape is correct.

x0 = 1;
xT = 0.5;
T = 6;
dt = 0.1;

t = 0:dt:T;
n = length(t);

x_true = zeros(1,61);
for i=1:n
    C1 = (0.5 + cos(24)/3 - 4*cos(6)/3)/sin(6);
    C2 = 4/3;
    x_true(i) = C1*sin(t(i)) + C2*cos(t(i)) - cos(4*t(i))/3;
end

v = -2*ones(n-2, 1);
u = ones((n-3), 1);
A = (diag(v) + diag(u,1) + diag(u,-1))/dt^2 + eye(n-2);

A9 = A;

b=zeros(n-2,1);
for i=1:n-2
    b(i) = 5*cos(4*t(i+1));
end
b(1) = b(1) - x0/dt^2;
b(end) = b(end) - xT/dt^2;

A10 = b;

x_int = A\b;
x = [x0;x_int;xT];

%plot(t,x_true,t,x,'ko','markersize',5)

A11 = x;
A12 = max(abs(x_true - A11'));

%%% Part 2: Shooting Method via Bisection
%%% Use the shooting method to solve the BVP
%%% It's a good idea to test out a few in the command window first to make
%%% sure that your initial conditions gets you to different sides of the right
%%% boundary condition.
%%% Use the plot to help you figure out what your choices of initial
%%% conditions should be

x_true = zeros(1,61);
for i=1:n
    C1 = (0.5 + cos(24)/3 - 4*cos(6)/3)/sin(6);
    C2 = 4/3;
    x_true(i) = C1*sin(t(i)) + C2*cos(t(i)) - cos(4*t(i))/3;
end

x0 = 1;
v1 = 2;
v2 = 4;
vmid = (v1+v2)/2;

[T, Y] = ode45(@(t, x) ddx(t,x),t,[x0,v1]);
xa = Y(:, 1);
[T, Y] = ode45(@(t, x) ddx(t,x),t,[x0,v2]);
xb = Y(:, 1);
[T, Y] = ode45(@(t, x) ddx(t,x),t,[x0,vmid]);
xmid = Y(:, 1);

%[T,Y] = ode45(@(t,x) ddx(t,x),t,[x0,v1]);
%plot(t,x_true,T,Y(:,1))

for i=1:100
    if abs(xmid(end) - xT) < 10e-8
        break
    else if sign(xmid(end) - xT) == sign(xa(end) - xT)
        v1 = vmid;
        [T,Y] = ode45(@(t,x) ddx(t,x),t,[x0,v1]);
        ta = T;
        xa = Y(:,1);
    else
        v2 = vmid;
        [T,Y] = ode45(@(t,x) ddx(t,x),t,[x0,v2]);
        tb = T;
        xb=Y(:,1);
        end
    end
    vmid=(v1+v2)/2;
    [T, Y] = ode45(@(t, x) ddx(t,x),t,[x0,vmid]);
    tmid = T;
    xmid = Y(:,1);
    %plot(t,x_true,tmid,xmid,'ko','markersize',3)
end

A13 = xmid;
A14 = max(abs(xmid' - x_true));
A15 = max(abs(A13-A11));

%%% You can set up your ODEs as functions here if you like

function dv = FHN(t,v, params)
    alpha = params(1);
    beta = params(2);
    tau = params(3);   
        
    dv1 = v(1) - ((v(1))^3)/3 - v(2) + (1/10)*(5+sin(pi*t/10));
    dv2 = (alpha + v(1) - beta*v(2))/tau;  
           
    dv = [dv1;dv2];
end

function dx = Chuas(t,x,params)
    alpha = params(1);
    beta = params(2);

    dx1 = alpha*(x(2) + x(1)/6 - ((x(1))^3)/16);
    dx2 = x(1) - x(2) + x(3);
    dx3 = -1*beta*x(2);
    
    dx = [dx1;dx2;dx3];
end

function dx = ddx(t,x)
    dx1 = x(2);
    dx2 = 5*cos(4*t) - x(1);

    dx = [dx1;dx2];
end