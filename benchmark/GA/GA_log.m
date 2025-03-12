function B = GA_log(M)
% B = -2log(M)
% M 是李群元素
% B 是恒等元切空间中对应的李代数 bivector
% 20M 5A 1D 1acos 1sqrt
if M(1) == 1
    B = [M(2:4); zeros(3, 1)];
else
    a = 1/(1 - M(1)^2);      % inv squared length
    b = acos(M(1))*sqrt(a);  % rotation scale
    c = a*M(8)*(1 - M(1)*b); % translation scale
    B = [c*M(5:7) + b*M(2:4); b*M(5:7)];
end
B = -2*B;
end