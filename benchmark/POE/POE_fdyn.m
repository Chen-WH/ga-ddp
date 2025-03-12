function ddq = POE_fdyn(rbt, q, dq, tau, fext)
if nargin == 4
    fext = zeros(6, 1);
end
n = rbt.n;
T0 = rbt.Tj;
L = rbt.Lj;
I = rbt.Ij;
type = rbt.type;

% init
invT = zeros(4, 4, n);
V = zeros(6, n);
c = zeros(6, n);
b = zeros(6, n);
d = zeros(1, n);
hI = zeros(6, 6, n);
hb = zeros(6, n);
dV = zeros(6, n);
ddq = zeros(n, 1);

% Forward iterations 速度与偏置项
invT(:, :, 1) = SE3_inv(T0(:, :, 1)*POE_motor(L(:, 1), q(1), type(1)));
V(:, 1) = dq(1)*L(:, 1);
c(:, 1) = SE3_adxi(V(:, 1))*L(:, 1)*dq(1);
b(:, 1) = -SE3_adxi(V(:, 1)).'*I(:, :, 1)*V(:, 1);
for i = 2:n
    invT(:, :, i) = SE3_inv(T0(:, :, i)*POE_motor(L(:, i), q(i), type(i)));
    V(:, i) = SE3_AdT(invT(:, :, i))*V(:, i-1) + dq(i)*L(:, i);
    c(:, i) = SE3_adxi(V(:, i))*L(:, i)*dq(i);
    b(:, i) = -SE3_adxi(V(:, i)).'*I(:, :, i)*V(:, i);
end

% Backward iterations 铰接体惯性与偏置力
for i = n:-1:2
    d(i) = 1/(L(:, i).'*I(:, :, i)*L(:, i));
    hI(:, :, i) = I(:, :, i) - d(i)*I(:, :, i)*L(:, i)*L(:, i).'*I(:, :, i);
    hb(:, i) = b(:, i) + hI(:, :, i)*c(:, i) + d(i)*I(:, :, i)*L(:, i)*(tau(i) - L(:, i).'*b(:, i));
    I(:, :, i-1) = I(:, :, i-1) + SE3_AdT(invT(:, :, i)).'*hI(:, :, i)*SE3_AdT(invT(:, :, i));
    b(:, i-1) = b(:, i-1) + SE3_AdT(invT(:, :, i)).'*hb(:, i);
end
d(1) = 1/(L(:, 1).'*I(:, :, 1)*L(:, 1));

% Forward iterations 加速度
dV(:, 1) = SE3_AdT(invT(:, :, 1))*[0; 0; 0; -rbt.g] + c(:, 1);
ddq(1) = d(1)*(tau(1) - L(:, 1).'*b(:, 1) - L(:, 1).'*I(:, :, 1).'*dV(:, 1));
dV(:, 1) = dV(:, 1) + L(:, 1)*ddq(1);
for i = 2:n
    dV(:, i) = SE3_AdT(invT(:, :, i))*dV(:, i-1) + c(:, i);
    ddq(i) = d(i)*(tau(i) - L(:, i).'*b(:, i) - L(:, i).'*I(:, :, i).'*dV(:, i));
    dV(:, i) = dV(:, i) + L(:, i)*ddq(i);
end
end