function out = grad_func(f,x,dx)
    out = zeros(size(x));
    s = size(x,1);
    for i=1:s
        A = zeros([s,1]);
        A(i,1) = 1;
        out(i,1) = (f(x + A*dx) - f(x - A*dx))/(2*dx);
    end
end