function ddq = GA_fdyn(rbt, q, dq, tau, fext)
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
c = zeros(6, n);
b = zeros(6, n);
d = zeros(1, n);
hI = zeros(6, 6, n);
hb = zeros(6, n);
dV = zeros(6, n);
ddq = zeros(n, 1);

% Forward iterations 速度与偏置项
M(:, 1) = GA_prodM(M0(:, 1), GA_motor(L(:, 1), q(1), type(1)));
Mt(:, 1) = GA_Rev(M(:, 1));
V(:, 1) = dq(1)*L(:, 1);
b(:, 1) = -PGA_Lie2(I(:, :, 1)*V(:, 1), V(:, 1));
for i = 2:n
    M(:, i) = GA_prodM(M0(:, i), GA_motor(L(:, i), q(i), type(i)));
    Mt(:, i) = GA_Rev(M(:, i));
    V(:, i) = PGA_Frame2(V(:, i-1), Mt(:, i)) + dq(i)*L(:, i);
    c(:, i) = dq(i)*PGA_Lie2(V(:, i), L(:, i));
    b(:, i) = -PGA_Lie2(I(:, :, i)*V(:, i), V(:, i));
end

% Backward iterations 铰接体惯性与偏置力
for i = n:-1:2
    d(i) = 1/(GA_metric(L(:, i))*I(:, :, i)*L(:, i));
    hI(:, :, i) = I(:, :, i) - d(i)*I(:, :, i)*L(:, i)*GA_metric(L(:, i))*I(:, :, i);
    hb(:, i) = b(:, i) + hI(:, :, i)*c(:, i) + d(i)*I(:, :, i)*L(:, i)*(tau(i) - GA_metric(L(:, i))*b(:, i));
    I(:, :, i-1) = I(:, :, i-1) + GA_AdM(M(:, i))*hI(:, :, i)*GA_AdM(Mt(:, i));
    b(:, i-1) = b(:, i-1) + PGA_Frame2(hb(:, i), M(:, i));
end
d(1) = 1/(GA_metric(L(:, 1))*I(:, :, 1)*L(:, 1));

% Forward iterations 加速度
dV(:, 1) = PGA_Frame2([-rbt.g; 0; 0; 0], Mt(:, 1)) + c(:, 1);
ddq(1) = d(1)*(tau(1) - GA_metric(L(:, 1))*(b(:, 1) + I(:, :, 1)*dV(:, 1)));
dV(:, 1) = dV(:, 1) + ddq(1)*L(:, 1);
for i = 2:n
    dV(:, i) = PGA_Frame2(dV(:, i-1), Mt(:, i)) + c(:, i);
    ddq(i) = d(i)*(tau(i) - GA_metric(L(:, i))*(b(:, i) + I(:, :, i)*dV(:, i)));
    dV(:, i) = dV(:, i) + ddq(i)*L(:, i);
end
end