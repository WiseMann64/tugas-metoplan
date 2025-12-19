function [x0,iter] = cg_polak_ribiere(f,x0)
    %% Conjugate Gradient
    g = @(x) grad_func(f,x,1e-6);
    c1 = 0.0001;
    c2 = 0.1;
    max_iter = 1000;
    
    tic
    eps = 1e-8;
    r = g(x0);
    p = -r;
    
    for k=1:max_iter
        if norm(r) < eps
            break
        end

        a0 = 0;
        amax = 1;
        a1 = 0.5;
        a = 0.5;
        phi = @(a) f(x0 + a*p);
        phip = @(a) grad_func(phi,a,1e-6);
    
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
    
        x = x0 + a*p;
        r = g(x);
        b = (g(x)'*(g(x) - g(x0)))/(norm(g(x0))^2);
        b = max([0 b]);
        p = -r + b*p;
    
        x0 = x;
    end
    toc
    fprintf("Found local minimum at x = [%s] with minimum value of %.6f\n",strjoin(compose("%.4f",x0), ";"),f(x0))
    iter = k;
end