
%%%  Problem 1
%%%  First go to the end of your m-file to create your Jacobi and
%%%  Gauss-Seidel functions.

%%% Once you have created your functions return here.
%%% Initialize your matrix A and RHS b
    
A=[1.1,0.2,-0.2,0.5;0.2,0.9,0.5,0.3;0.1,0,1,0.4;0.1,0.1,0.1,1.2];
b=[1;0;1;0];

%%% Use your Jacobi and Gauss-Seidel functions to find A1 through A4.
[T,E]=jacobi(A,b,1e-2);
Tj_2=T;
Ej_2=max(abs(E));
[T,E]=jacobi(A,b,1e-4);
Tj_4=T;
Ej_4=max(abs(E));
[T,E]=jacobi(A,b,1e-6);
Tj_6=T;
Ej_6=max(abs(E));
[T,E]=jacobi(A,b,1e-8);
Tj_8=T;
Ej_8=max(abs(E));

A1=[Tj_2, Tj_4, Tj_6, Tj_8];
A2=[Ej_2, Ej_4, Ej_6, Ej_8];

[T,E]=gauss(A,b,1e-2);
Tgs_2=T;
Egs_2=max(abs(E));
[T,E]=gauss(A,b,1e-4);
Tgs_4=T;
Egs_4=max(abs(E));
[T,E]=gauss(A,b,1e-6);
Tgs_6=T;
Egs_6=max(abs(E));
[T,E]=gauss(A,b,1e-8);
Tgs_8=T;
Egs_8=max(abs(E));
A3=[Tgs_2, Tgs_4, Tgs_6, Tgs_8];
A4=[Egs_2, Egs_4, Egs_6, Egs_8];

%%%  Problem 2
%%%  Initialize your Day 0 vector x
    
x=[0.9;0.09;0.01];

      
%%%  Part 1: without a vaccine
%%%  Make sure to have p = 0
%%%  Initialize the SIR matrix M, and save it as A5
p=0;
M=[1,0,0;0,1,0;0,0,1] + [-((1/200)+p),0,1/10000;1/200,-1/1000,0;p,1/1000,-1/10000];
A5=M;

%%%  Create a loop to find the day that the number of infected
%%%  individuals hits 50% and another loop for the steady state of the
%%%  infected population
%%%  There is a way to put everything under one loop if you make clever use
%%%  of conditionals

day=0;
for n=1:10000
    x=M*x;
    day=day+1;
    if x(2)>=0.5
        break
    end 
end
D0=day;

x=[0.9;0.09;0.01];
day=0;
for n=1:100000
    x=M*x;
    day=day+1;
end
F0=x(2);

%%% Save the days and steady state in a row vector A6

A6=[D0,F0];

%%%  Reinitialize your Day 0 vector x

x=[0.9;0.09;0.01];
     

%%%  Part 2: with a vaccine
%%%  Make sure to have p = 2/1000
%%%  Initialize the SIR matrix M, and save it as A7

p=2/1000;
M=[1,0,0;0,1,0;0,0,1] + [-((1/200)+p),0,1/10000;1/200,-1/1000,0;p,1/1000,-1/10000];
A7=M;


%%%  Create a loop to find the day that the number of infected
%%%  individuals hits 50% and another loop for the steady state of the
%%%  infected population
%%%  There is a way to put everything under one loop if you make clever use
%%%  of conditionals

day=0;
for n=1:10000
    x=M*x;
    day=day+1;
    if x(2)>=0.5
        break
    end 
end
D1=day;

x=[0.9;0.09;0.01];
day=0;
for n=1:100000
    x=M*x;
    day=day+1;
end
F1=x(2);


%%% Save the days and steady state in a row vector A8

A8=[D1,F1];
 
 
%%%  Problem 3
  
%%%  Initialize your 114x114 tridiagonal matrix A

A=zeros(114,114);
A=A+diag(2-diag(A));
A=A+diag((-1-diag(A,-1)),-1);
A=A+diag((-1-diag(A,+1)),+1);

A9=A;

%%%  Initialize your 114x1 RHS column vector rho

rho=zeros(114,1);
for j=1:114
    rho(j)=2*(1-(cos(53*pi/115)))*(sin(53*pi*j/115));
end

A10=rho;

%%%  Implement Jacobi's method for this system.
%%%  Don't use the function you created before because that was designed for
%%%  a specific task, and will not work here.

L=tril(A,-1);
U=triu(A,+1);
D=diag(A);
M=-(L+U)./D;
c=rho./D;
phi_old=ones(114,1);
iter=1;
for n=1:1e5
    phi_new=M*phi_old + c;
    E=max(abs(phi_new-phi_old));
    phi_old=phi_new;
    iter=iter+1;
    if E<1e-5
        break
    end
end

A11=phi_new;
A12=iter;

%%%  Create a column vector phi that contains the exact solution given in
%%%  the assignment file

phi=zeros(114,1);
for j=1:114
    phi(j)=(sin(53*pi*j/115));
end

%%%  Save the difference of the Jacobi solution and the exact solution as
%%%  A13.  Use the maximal entry in absolute value to calculate this error.

A13=max(abs(phi-A11));

%%%  Implement Gauss-Seidel for this system.
%%%  Don't use the function you created before because that was designed for
%%%  a specific task, and will not work here.

LpD=tril(A);
U=triu(A,+1);
M=-LpD\U;
c=LpD\rho;
phi_old=ones(114,1);
iter=1;
for n=1:1e5
    phi_new=M*phi_old + c;
    E=max(abs(phi_new-phi_old));
    phi_old=phi_new;
    iter=iter+1;
    if E<1e-5
        break
    end
end

A14=phi_new;
A15=iter;

%%%  Save the difference of the Gauss-Seidel solution and the exact solution as
%%%  A13.  Use the maximal entry in absolute value to calculate this error.


A16=max(abs(phi-A14));


%%% Jacobi and Gauss Seidel Iteration functions
%%% Create your functions here
%%% Both functions will need two outputs and three inputs
%%% The code within the function will be very similar to
%%% Week 4 coding lecture 2


function [T,E] = jacobi(A,b,eps)
    x = zeros(4,1); 
    L = tril(A,-1);
    U = triu(A,+1);
    D = diag(A);
    M = -(L+U)./D;
    c = b./D;
    N = 1000;
    T = 0;
    for n=1:N
        x=M*x + c;
        E=A*x - b;
        T=T+1;
        if (max(abs(E))) < eps
            break
        end
    end    
end

function [T,E] = gauss(A,b,eps)
    x = zeros(4,1); 
    LpD = tril(A);
    U = triu(A,+1);
    M = -LpD\U;
    c = LpD\b;
    N = 1000;
    T = 0;
    for n=1:N
        x=M*x + c;
        E=A*x - b;
        T=T+1;
        if (max(abs(E))) < eps
            break
        end
    end    
end