
A = readmatrix('bridge_matrix.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Don't delete anything above this line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Problem 1
%%% A is already initialized above from a separate file (don't delete line 4).
%%% Initialize the data (right hand side) b.

b=[0;0;0;0;0;0;0;0;3;0;(exp(1)^2);0;pi];

%%% Solve for the force vector and save it as A1

x=A\b;
A1=x;

%%% Compute the PA = LU factorization of A
%%% Remember MATLAB computes the PA = LU factorization as A = PLU.
%%% Save L AS A2, and c as A3.

[L,U,P]=lu(A);
A2=L;
c=L\(P*b);
A3=c;

%%% Create a loop that breaks when one of the forces is greater than 20 tons
%%% Save A4 as the weight of the truck in position 8
%%% Save A5 as the entry of the force vector that exceeds 20 tons

numb=1e5;
for n=1:numb
    b(9)=b(9)+0.001;
    x=A\b;
    if max(abs(x))>=20
        break
    end
end
A4=b(9);
A5=find(abs(x)>=20);

%%% Problem 2
%%% Initialize, alpha, omega, and A, and compute the PA = LU factorization
%%% Remember MATLAB computes the PA = LU factorization as A = PLU.

alpha=-0.002;
omega=0.06;
A=[1-alpha,-omega;omega,1-alpha];
[L,U,P]=lu(A);


%%% The initializations can get a little tricky so definitely ask for help
%%% if you're stuck.
%%% Initialize a matrix made up of the position vector at each time
%%% Set the first x and y coordinates at time = 0 in your matrix
%%% to the values instructed in the assignment file.
%%% Create a loop that loops through each time given in the assignment file.
%%% Compute the new right hand side c using P, L, and/or U.
%%% You may need to recall that the inverse of P is P transpose
%%% Solve for the position by solving the Ux = c equation.
%%% Save all x coordinates as A6
%%% Save all y coordinates as A7
%%% Save the distance from the origin as A8

x=zeros(2,1001);
x(:,1)=[1;-1];
for n=2:1001
    x(:,n)=A\x(:,n-1);
end
A6=x(1,:);
A7=x(2,:);
A8=(A6.^2 + A7.^2).^(1/2);


%%% Initialize a position vector
%%% Initialize a distance variable
%%% Initialize a time variable
%%% Create a loop that breaks when the distance from the origin is
%%% less than 0.06.
%%% In the loop compute the position using P, L, and/or U and
%%% compute the distance from the origin.
%%% Iterate time at each iteration of the loop.
%%% Save the time the loop breaks as A9.
%%% Save the distance from the origin as A10.
numb=1001;
x=zeros(2,numb);
x(:,1)=[1;-1];
for n=2:numb
    x(:,n)=A\x(:,n-1);
    dist=((x(1,:)).^2 + (x(2,:)).^2).^(1/2);
    if dist(n)<0.06
        break
    end
end
A9=find(dist<0.06,1,'first')-1;
A10=dist(find(dist<0.06,1,'first'));

%%% Problem 3
%%% First go to the end of the file to create your function

%%% After you make your function come back to this line.
%%% Save A11 as R(pi/8)
%%% Rotate the vector given in the assignment file and save it as A12.
theta=pi/8;
A11=rotation_R(pi/8);
x=[pi/10;2.1;-exp(1)];
R=rotation_R(pi/3);
A12=R*x;

%%% Find the vector x that was rotated to give you vector b.
%%% Save the vector x as A13

b=[1.4;-pi/10;2.8];
R=rotation_R(pi/6);
A13=R\b;


%%% Invert the R(3*pi/4) and save it as A14.
%%% Find the angle theta that would give you this inverse
%%% without having to do matrix operations, and save the angle
%%% as A15.

R=rotation_R(3*pi/4);
A14=inv(R);
A15=-3*pi/4;



%%% Create a function here for the rotation matrix that
%%% takes an input in radians and outputs the matrix.
%%% After you make the function go back to Problem 3

function R = rotation_R(theta)
    for R=zeros(3,3)
        R=[(cos(theta)),0,(sin(theta));0,1,0;(-sin(theta)),0,(cos(theta))];
    end
end
