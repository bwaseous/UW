%%% Problem 1
%%% Initialize t, and x_true

dt = 0.1;
t = (0:dt:10);
x_true = (1/2)*(cos(t) + sin(t) + exp(-t));

%%% Forward Euler
%%% Write a forward Euler scheme

N = length(t);
x = zeros(1,N);
x0 = 1;
x(1) = x0;

for n = 1:N-1
    x(n+1) = x(n) + dt*(cos(t(n)) - x(n));
end

%plot(t,x,'k',t,x_true)

A1 = x;

A2 = abs(x - x_true);

%%% Backward Euler
%%% Write a backward Euler scheme

N = length(t);
x = zeros(1,N);
x0 = 1;
x(1) = x0;

for n = 1:N-1
    x(n+1) = (x(n) + dt*cos(t(n+1))) / (1 + dt);
end

A3 = x;

A4 = abs(x - x_true);

%%% Built-in Solver
%%% Use ode45
%%% Don't forget to transpose the solution you get from ode45.

tspan = (0:0.1:10);
x0 = 1;
[t, x] = ode45(@(t,x) myODE(t, x) ,tspan, x0);

A5 = x';

x_true = (1/2)*(cos(t) + sin(t) + exp(-t));


A6 = abs(x' - x_true');

%%% Problem 2
%%% Initialize the parameters

a = 8;

%%% Forward Euler for dt = 0.01
%%% Initialize t and x_true
%%% Write a forward Euler scheme

dt = 0.01;
t = (0:dt:2);
x_true = 2*atan((exp(a*t)) / (1 + sqrt(2)));
N = length(t);
x = zeros(1,N);
x0 = pi/4;
x(1) = x0;

for n = 1:N-1
        x(n+1) = x(n) + dt*(a*sin(x(n)));
end

%plot(t,x_true, t, x)

A7 = x;

A8 = max(abs(x - x_true));

%%% Forward Euler dt = 0.001
%%% Reinitialize t and x_true
%%% Write a forward Euler scheme

dt = 0.001;
t = (0:dt:2);
x_true = 2*atan((exp(a*t)) / (1 + sqrt(2)));
N = length(t);
x = zeros(1,N);
x0 = pi/4;
x(1) = x0;

for n = 1:N-1
        x(n+1) = x(n) + dt*(a*sin(x(n)));
end

%plot(t,x_true, t, x)

A9 = A8 / max(abs(x - x_true));

%%% Predictor-Corrector dt = 0.01
%%% Reinitialize t and x_true
%%% Write a forward Euler scheme and a backward Euler scheme in the same
%%% loop.
%%% The forward Euler scheme is the predictor.  The answer from forward
%%% Euler will be plugged into the sin(x_n+1) in the backward Euler scheme.

a=8;
dt = 0.01;
t = (0:dt:2);
x_true = 2*atan((exp(a*t)) / (1 + sqrt(2)));
N = length(t);
x = zeros(1,N);
x0 = pi/4;
x(1) = x0;
xtilda = x0;

for n = 1:N-1
        xtilda = x(n) + dt*(a*sin(xtilda));
        x(n+1) = x(n) + dt*(a*sin(xtilda));
        xtilda = x(n+1);
end

%plot(t,x_true, t, x)

A10 = x;

A11 = max(abs(x - x_true));

%%% Predictor-Corrector dt = 0.001
%%% Reinitialize t and x_true
%%% Write a forward Euler scheme and a backward Euler scheme in the same
%%% loop.
%%% The forward Euler scheme is the predictor.  The answer from forward
%%% Euler will be plugged into the sin(x_n+1) in the backward Euler scheme.

a=8;
dt = 0.001;
t = (0:dt:2);
x_true = 2*atan((exp(a*t)) / (1 + sqrt(2)));
N = length(t);
x = zeros(1,N);
x0 = pi/4;
x(1) = x0;
xtilda = x0;

for n = 1:N-1
        xtilda = x(n) + dt*(a*sin(xtilda));
        x(n+1) = x(n) + dt*(a*sin(xtilda));
        xtilda = x(n+1);
end

%plot(t,x_true, t, x)

A12 = A11 / max(abs(x - x_true));

%%% Builtin Solver
%%% Reinitialize t and x_true
%%% Use ode45 to solve the ODE.
%%% Don't forget to transpose the solution you get from ode45.

tspan = (0:0.01:2);
x0 = pi/4;
a = 8;

[t, x] = ode45(@(t,x) myODE2(t, x) ,tspan, x0);

A13 = x';

x_true = 2*atan((exp(a*t)) / (1 + sqrt(2)));

A14 = max(abs(x' - x_true'));

tspan = (0:0.001:2);
x0 = pi/4;
a = 8;

[t, x] = ode45(@(t,x) myODE2(t, x) ,tspan, x0);

x_true = 2*atan((exp(a*t)) / (1 + sqrt(2)));

A15 = A14 / max(abs(x' - x_true'));

%%%  If you want to write local functions, put them here

function xn = myODE(t, x0)
    xn = cos(t) - x0;
end

function xn2 = myODE2(t, x0)
    a = 8;
    xn2 = a*sin(x0);
end