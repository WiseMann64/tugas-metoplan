function [x0,iter] = dogleg(f,x0)
    h = 1e-6;
    fx = @(x) (f(x+[h;0]) - f(x-[h;0]))/(2*h);
    fy = @(x) (f(x+[0;h]) - f(x-[0;h]))/(2*h);
    
    fxx = @(x) (fx(x+[h;0]) - fx(x-[h;0]))/(2*h);
    fxy = @(x) (fx(x+[0;h]) - fx(x-[0;h]))/(2*h);
    fyy = @(x) (fy(x+[0;h]) - fy(x-[0;h]))/(2*h);
    
    g = @(x) [fx(x);fy(x)];
    H = @(x) [fxx(x),fxy(x);fxy(x),fyy(x)];
    
    max_iter = 1000;
    Dm = 1;
    D0 = 0.5;
    gamma = 0.1;
    
    tic
    for k=1:max_iter
        H0 = H(x0);
        H0 = (H0 + H0')/2;
        g0 = g(x0);
    
        if norm(g0) < 1e-8
            break
        end
    
        hsd = -g0;
        hgn = -H0\g0;
        alpha = (g0'*g0)/(g0'*H(x0)*g0);
    
        if norm(hgn) <= D0
            p = hgn;
        elseif norm(alpha*hsd) >= D0
            p = D0/norm(alpha*hsd)*hsd;
        else
            func = @(b) norm(alpha*hsd + b*(hgn - alpha*hsd)) - D0;
            options = optimset("Display","off");
            b = fsolve(func,1,options);
            p = alpha*hsd + b*(hgn - alpha*hsd);
        end
        
        pred = -g0'*p - 0.5*p'*H0*p;
        ared = f(x0) - f(x0 + p);
        r = ared/pred;
        if r < 1/4
            D0 = norm(p)/4;
        elseif r > 3/4 && abs(norm(p) - D0) < 1e-6
            D0 = min([2*D0 Dm]);
        end
        if r > gamma
            x0 = x0 + p;
        end
    end
    toc
    fprintf("Found local minimum at x = [%s] with minimum value of %.6f\n",strjoin(compose("%.4f",x0), ";"),f(x0))

    iter = k;
end