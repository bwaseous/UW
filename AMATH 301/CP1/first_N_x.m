function x = first_N_x(r,x0,N)
    x(1)=x0;
    for n=2:N
        x(n)=r*(x(n-1))*(1-(x(n-1)));
    end
end

