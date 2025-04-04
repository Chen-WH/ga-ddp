function J = PGA_jacob3(rbt, q, Q0)
n = rbt.n;
M0 = rbt.Mj;
L = rbt.Lj;
type = rbt.type;
J = zeros(4, n);
Mf = zeros(8, n);
Mb = zeros(8, n);
Mf(:, 1) = GA_prodM(M0(:, 1), GA_motor(L(:, 1), q(1), type(1)));
for i = 2:n
    Mf(:, i) = GA_prodM(Mf(:, i-1), M0(:, i), GA_motor(L(:, i), q(i), type(i)));
end
Mb(:, n) = M0(:, end);
for i = n:-1:2
    Mb(:, i-1) = GA_prodM(M0(:, i), GA_motor(L(:, i), q(i), type(i)), Mb(:, i));
end
for i = 1:n
    J(:, i) = PGA_Frame3(PGA_Lie3(PGA_Frame3(Q0, Mb(:, i)), L(:, i)), Mf(:, i));
end
end