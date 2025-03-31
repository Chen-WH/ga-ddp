function ad_B = GA_adB(B)
% Optimized for vectorized computation
% ad_A(B) = [A, B] = -(AB - BA)/2
% basis vector = ["1", "e0", "e1", "e2", "e3", "e01", "e02", "e03", "e23", "e31", "e12", "e123", "e023", "e031", "e012", "e0123"]
n = size(B, 2);
if n == 1
    ad_B = [0, -B(6),  B(5),     0, -B(3),  B(2);
         B(6),     0, -B(4),  B(3),     0, -B(1);
        -B(5),  B(4),     0, -B(2),  B(1),     0;
            0,     0,     0,     0, -B(6),  B(5);
            0,     0,     0,  B(6),     0, -B(4);
            0,     0,     0, -B(5),  B(4),    0];
else
    ad_B = zeros(6, 6, n);
    for num = 1:n
        ad_B(:, :, num) = [0, -B(6, num),  B(5, num),          0, -B(3, num),  B(2, num);
                   B(6, num),          0, -B(4, num),  B(3, num),          0, -B(1, num);
                  -B(5, num),  B(4, num),          0, -B(2, num),  B(1, num),          0;
                  0,          0,          0,          0, -B(6, num),  B(5, num);
                  0,          0,          0,  B(6, num),          0, -B(4, num);
                  0,          0,          0, -B(5, num),  B(4, num),         0];
    end
end
end