
%%% Problem 1
data = readmatrix('lynx.csv');
t = data(1, :);
pop = data(2, :);
%%% Don't delete the lines above when submitting to gradescope

%%% Replace the value of the population for the years given in the assignment file and save it as A1

pop(11) = 34;
pop(29) = 27;
A1= pop;

%%% Calculate the value of the cubic spline interpolation of the data at t = 24.5 using the interp1 function.  Save this as A2.
xpop = (0:0.001:length(pop));
ypop = interp1(t, pop, xpop, 'spline');
%plot(t , pop, 'ko', xpop, ypop)
A2 = ypop(24501);

%%% Use polyfit to calculate the coefficients for A3, A5, and A7
%%% Use norm to calculate the error for A4, A6, and A8

coeffs = polyfit(t, pop, 1);
A3 = coeffs;
y = polyval(coeffs,t);
A4 = norm(y-pop);

coeffs = polyfit(t, pop, 2);
A5 = coeffs;
y = polyval(coeffs,t);
A6 = norm(y-pop);

coeffs = polyfit(t, pop, 10);
A7 = coeffs;
y = polyval(coeffs,t);
A8 = norm(y-pop);

%%% Problem 2
data = readmatrix('CO2_data.csv');
t = data(1, :);
co2 = data(2, :);
%%% Don't delete the lines above when submitting to gradescope

%%% Use polyfit to calculate the coefficients for A9
%%% Use norm to calculate the error for A10

coeffs = polyfit(t, co2, 1);
A9 = coeffs;
f = polyval(coeffs,t);
A10 = norm(f-co2);

%%% Fit the exponential

co2_b = co2-260;
coeffs = polyfit(t, log(co2_b), 1);
A11 = zeros(3,1); 
A11 = coeffs; 
A11(1)=exp(1)^coeffs(2); 
A11(2) = coeffs(1); 
A11(3) = 260;
f = (A11(1)*exp(1).^(A11(2).*t)) + 260;
A12 = norm(f - co2);

%%% Fit the sinusoidal
%%% There are a few different ways to do this, and we will refrain from giving away the answer to this part.  The class has been doing loops for a while now, so this part should be doable, albeit a little tricky.  We can however check to see if there are any bugs that we can spot.

co2_744 = co2(1:744);
t_744 = t(1:744);
co2_sin = co2_744 - ((A11(1)*exp(1).^(A11(2).*t_744)) + 260);

Amp = zeros(1,62);
month = 0;
m = 0;
for i=1:62
    month = co2_sin(1,(12*m + 1):(12+(12*m)));
    min_m = min(month);
    max_m = max(month);
    Amp(i) = (max_m - min_m)./2;
    m = m+1;
end

Amp_avg = (sum(Amp)/62);
B = 2*pi;

A13 = [Amp_avg, B];
f = ((A11(1)*exp(1).^(A11(2).*t)) + 260) + (Amp_avg.*sin((B).*t));
A14 = norm(f-co2);

plot(t,f,t,co2,'ko','markersize',2)
%plot(t,co2_sin)