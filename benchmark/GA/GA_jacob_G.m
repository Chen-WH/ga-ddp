function J = GA_jacob_G(rbt, q)
n = rbt.n;
M0 = rbt.Mj;
L = rbt.Lj;
type = rbt.type;
J = zeros(6, n);
M = M0(:, 1);
J(:, 1) = -PGA_Frame2(L(:, 1), M)/2;
for i = 1:rbt.n-1
    M = GA_prodM(M, GA_motor(L(:, i), q(i), type(i)), M0(:, i+1));
    J(:, i+1) = -PGA_Frame2(L(:, i+1), M)/2;
end
end