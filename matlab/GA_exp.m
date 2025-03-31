function M = GA_exp(B)
% M = exp(-B/2)
% B 是李代数 bivector
% M 是 Motor Group 上对应的李群元素
% 23M 8A 2D 1sin 1cos 1sqrt

B = -0.5*B;
l = B(4:6)'*B(4:6);
if l == 0
    M = [1; B(1:3); 0; 0; 0; 0];
else
    m = B(1:3)'*B(4:6);
    a = sqrt(l);
    c = cos(a);
    s = sin(a)/a;
    t = m/l*(c - s);
    M = [c; s*B(1:3) + t*B(4:6); s*B(4:6); m*s];
end
end