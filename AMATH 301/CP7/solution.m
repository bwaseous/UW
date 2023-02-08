%%% Problem 1
data = readmatrix('population.csv');
t = data(1, :);
N = data(2, :);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Determine your stepsize dt from the vector t

dt = 1;
deltat = 10;
%plot(t,N,'ko--','markersize',3)

%%% Use the appropriate second order differences from the Theory Lecture

t_deriv = 24;
A1 = (((3.*N(t_deriv)) - (4.*N(t_deriv - dt)) + (N(t_deriv - 2*dt))))/(2*deltat);

t_deriv = 10;
A2 = ((N(t_deriv + dt)) - (N(t_deriv - dt)))/(2*deltat);

t_deriv = 1;
A3 = ((-3*N(t_deriv)) + (4*N(t_deriv + dt)) - (N(t_deriv + 2*dt)))/(2*deltat);

%%% For dN/dt you will need to use a combination of the above differences,
%%% but the choice will be obvious based on which direction you can/cannot
%%% go in the horizontal axis.  Whenever possible use central difference;
%%% only use forward or backward when central is not possible.

dN = zeros(1,24);
per_capita = zeros(1,24);
dt = 1;
deltat = 10;
for k=1:length(t)
    t_deriv = k;
    if k ~= 1 && k~= 24
        dN(k) = ((N(t_deriv + dt)) - (N(t_deriv - dt)))/(2*deltat);
        per_capita(k) = dN(k)./N(k);
    elseif k == 1
        dN(k) = ((-3*N(t_deriv)) + (4*N(t_deriv + dt)) - (N(t_deriv + 2*dt)))/(2*deltat);
        per_capita(k) = dN(k)./N(k);
    elseif k == 24
        dN(k) = (((3.*N(t_deriv)) - (4.*N(t_deriv - dt)) + (N(t_deriv - 2*dt))))/(2*deltat);
        per_capita(k) = dN(k)./N(k);
    end
end

A4 = dN;
A5 = per_capita;
A6 = sum(per_capita)/length(per_capita);

%%% Problem 2
data = readmatrix('brake_pad.csv');
r = data(1, :);
T = data(2, :);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Determine your stepsize dr from the vector r

dr = 0.017; 
thetap = 0.7051;

%%% Use the LHR formula from the coding lecture

T_tot_LHR = 0;
for k  = 1:length(r)-1
    T_tot_LHR = T_tot_LHR + (r(k))*(T(k)*thetap*dr);
end

A7 = T_tot_LHR;

A_LHR = 0;
for k = 1:length(r)-1
    A_LHR = A_LHR + (r(k)*thetap*dr);
end

A8 = T_tot_LHR/A_LHR;

%%% Use the RHR formula from the coding lecture

T_tot_RHR = 0;
for k = 2:length(r)
    T_tot_RHR = T_tot_RHR + (r(k))*(T(k)*thetap*dr);
end

A9 = T_tot_RHR;

A_RHR = 0;
for k = 2:length(r)
    A_RHR = A_RHR + (r(k)*thetap*dr);
end

A10 = T_tot_RHR/A_RHR;


%%% Use the Trapezoid rule formula or the trapz function from the coding lecture

T_tot_trap = (T_tot_LHR + T_tot_RHR)/2;
A11 = T_tot_trap;

A_trap = (A_LHR + A_RHR)/2;
A12 = T_tot_trap/A_trap;


%%% Problem 3
%%% You'll have to use anonymous functions here.  You can see the syntax in
%%% the Numerical Integration coding lecture where the builtin function
%%% "integrand" is used.

Tf = @(x, mu) (mu ./ (sqrt((((mu).^2)./2 - ((mu).^3)./3) - (((mu.*x).^2)./2 - ((mu.*x).^3)./3))));
a=0; b=1;

A13 = integral(@(x) Tf(x,0.95),a,b);
A14 = integral(@(x) Tf(x,0.5),a,b);
A15 = integral(@(x) Tf(x,0.01),a,b);
