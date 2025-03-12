function ad = SE3_adxi(xi)
% given Lie algebra xi, calculate ad_xi
n = size(xi, 2);
if n == 1
    ad = [skew(xi(1:3)), zeros(3); skew(xi(4:6)), skew(xi(1:3))];
else
    ad = zeros(6, 6, n);
    for num = 1:n
        ad(:, :, num) = [skew(xi(1:3, num)), zeros(3); skew(xi(4:6, num)), skew(xi(1:3, num))];
    end
end
end