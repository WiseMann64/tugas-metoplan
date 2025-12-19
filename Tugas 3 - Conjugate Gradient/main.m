clc
% Griewank Function
f = @(x) sum(x.^2)/4000 - prod(cos(x./sqrt((1:size(x,1))'))) + 1;

[x,i] = cg_polak_ribiere(f,[4.3;8.1]);
f(x),i

[x,i] = cg_polak_ribiere(f,[0.2;0.3]);
f(x),i

[x,i] = cg_polak_ribiere(f,[4.5;9.0]);
f(x),i

[x,i] = cg_polak_ribiere(f,[2.3;5.7]);
f(x),i

[x,i] = cg_polak_ribiere(f,[3.2;6.8]);
f(x),i