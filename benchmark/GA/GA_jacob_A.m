function J = GA_jacob_A(rbt, q)
n = rbt.n;
M0 = rbt.Mj;
L = rbt.Lj;
type = rbt.type;
J = zeros(8, n);
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
    J(:, i) = GA_prodM(GA_prodMB(Mf(:, i), L(:, i)), Mb(:, i))/2;
end
end