function tau = GA_idyn(rbt, q, dq, ddq, fext)
if nargin == 4
    fext = zeros(6, 1);
end
n = rbt.n;
M0 = rbt.Mj;
L = rbt.Lj;
I = rbt.Ij;
type = rbt.type;

% init
M = zeros(8, n);
Mt = zeros(8, n);
V = zeros(6, n);
dV = zeros(6, n);
F = zeros(6, n);
tau = zeros(n, 1);

% Forward iterations
M(:, 1) = GA_prodM(M0(:, 1), GA_motor(L(:, 1), q(1), type(1)));
Mt(:, 1) = GA_Rev(M(:, 1));
V(:, 1) = dq(1)*L(:, 1);
dV(:, 1) = PGA_Frame2([-rbt.g; 0; 0; 0], Mt(:, 1)) + ddq(1)*L(:, 1);
for i = 2:n
    M(:, i) = GA_prodM(M0(:, i), GA_motor(L(:, i), q(i), type(i)));
    Mt(:, i) = GA_Rev(M(:, i));
    V(:, i) = PGA_Frame2(V(:, i-1), Mt(:, i)) + dq(i)*L(:, i);
    dV(:, i) = PGA_Frame2(dV(:, i-1), Mt(:, i)) + dq(i)*PGA_Lie2(V(:, i),L(:, i)) + ddq(i)*L(:, i);
end

% Backward iterations
F(:, n) = -PGA_Frame2(fext, M0(:, n+1)) + I(:, :, n)*dV(:, n) - PGA_Lie2(I(:, :, n)*V(:, n), V(:, n));
tau(n) = GA_metric(L(:, n))*F(:, n);
for i = n-1:-1:1
    F(:, i) = PGA_Frame2(F(:, i+1), M(:, i+1)) + I(:, :, i)*dV(:, i) - PGA_Lie2(I(:, :, i)*V(:, i), V(:, i));
    tau(i) = GA_metric(L(:, i))*F(:, i);
end
end