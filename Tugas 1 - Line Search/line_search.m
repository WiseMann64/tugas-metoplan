function [x0,iter] = line_search(x0)
%% Algoritma 3.5 dan 3.6
    % Griewank Function
    f = @(x) sum(x.^2)/4000 - prod(cos(x./sqrt((1:size(x,1))'))) + 1;
    
    % Gradient Function (numerical with central difference)
    g = @(x) grad_func(f,x,0.000001);
    
    % Parameter
    c1 = 0.0001;
    c2 = 0.1;
    max_iter = 1000;
    
    tic
    for k=1:max_iter
        % p menggunakan -grad F(x)
        p = -g(x0);
    
        % Line search
        a0 = 0;
        amax = 1;
        a1 = 0.5;
        a = 0.5;
        phi = @(a) f(x0 + a*p);
        phip = @(a) grad_func(phi,a,0.0001);
    
        for i=1:max_iter
            if phi(a1) > phi(0) + c1*a1*phip(0) || (i > 1 && phi(a1) >= phi(a0))
                a = zoom(a0,a1,phi,c1,c2);
                break
            end
            if abs(phip(a1)) <= -c2*phip(0)
                a = a1;
                break
            end
            if phip(a1) >= 0
                a = zoom(a1,a0,phi,c1,c2);
                break
            end
            a0 = a1;
            a1 = (a1 + amax)/2;
        end
        
        x0 = x0 + a*p;
        if norm(g(x0)) < 1e-8
            break
        end
    end
    toc
    fprintf("Found local minimum at x = [%s] with minimum value of %.6f\n",strjoin(compose("%.4f",x0), ";"),f(x0))
    fprintf("Process took %d iterations\n",k)
    iter = k;
end