%A=[1.1,0.2,-0.2,0.5;0.2,0.9,0.5,0.3;0.1,0,1,1.2;0.1,0.1,0.1,1.2]; b=[1;0;1;0];
function [T,E] = jacobi(A,b,eps)
    x = zeros(4,1); 
    L = tril(A,-1);
    U = triu(A,+1);
    D = diag(A);
    M = -(L+U)./D;
    c = b./D;
    N = 1000;
    T = 1;
    for n=1:N
        x=M*x + c;
        E=A*x - b;
        T=T+1;
        if (max(abs(E))) < eps;
            break
        end
    end    
end

