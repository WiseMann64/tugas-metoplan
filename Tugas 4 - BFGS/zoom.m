function a = zoom(a1,a2,phi,c1,c2)
    phip = @(x) grad_func(phi,x,0.0001);
    max_iter = 100;
    for i=1:max_iter
        % Bisection
        aj = (a1 + a2)/2;
        
        phiaj = phi(aj);
        
        if phiaj > phi(0) + c1*aj*phip(0) || phi(aj) >= phi(a1)
            a2 = aj;
        else
            if abs(phip(aj)) <= -c2*phip(0)
                a = aj;
                return
            end
            if phip(aj)*(a2 - a1) >= 0
                a2 = a1;
            end
            a1 = aj;
        end
    end
    a = a2;
end