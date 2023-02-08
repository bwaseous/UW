%% 7
A = [1,0,0;1,1,0;0,1,1];
b = [1;1.1;-0.9];
x = (A'*A)\(A'*b);
A1 = x;
A2 = sum(x);

P = (1/3).*[2, -1, -1; -1, 2, -1; -1, -1, 2];
xnaive = P*x;
A3 = xnaive;
A4 = (1/2)*(norm((A*xnaive) - b))^2;

z = ((A*P)'*A*P)\((A*P)'*b);
xcorr = P*z;
A5 = xcorr;
A6 = (1/2)*(norm((A*xcorr) - b))^2;


%% 8
A = [1, 0, 0, 1, 1; 0, 1, 1, 0, 0; 1, 0, 1, 0, 1; 1, 0, 0, 1, 1];
y = [1;2;3;4];
V = A';
v1 = V(:,1);
v2 = V(:,2);
v3 = V(:,3); 
%v4 = V(:,4);
q1 = v1/norm(v1);
v2t = v2 - dot(q1,v2)*q1;
q2 = v2t/norm(v2t);
v3t = v3 - q1*dot(q1,v3) - q2*dot(q2,v3);
q3 = v3t / norm(v3t);
%v4t = v4 - q1*dot(q1,v4) - q2*dot(q2,v4) - q3*dot(q3,v4);
%q4 = v4t / norm(v4t);
Q = [q1, q2, q3];
R = Q'*V;
xhat = Q*(R'\y);
A7 = xhat;
A8 = (1/2)*norm(A*xhat - y)^2;