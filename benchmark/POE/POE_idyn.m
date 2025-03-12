function tau = POE_idyn(rbt, q, dq, ddq, fext)
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
dV = zeros(6, n);
F = zeros(6, n);
tau = zeros(n, 1);

% Forward iterations
invT(:, :, 1) = SE3_inv(T0(:, :, 1)*POE_motor(L(:, 1), q(1), type(1)));
V(:, 1) = dq(1)*L(:, 1);
dV(:, 1) = SE3_AdT(invT(:, :, 1))*[0; 0; 0; -rbt.g] + ddq(1)*L(:, 1);
for i = 2:n
    invT(:, :, i) = SE3_inv(T0(:, :, i)*POE_motor(L(:, i), q(i), type(i)));
    V(:, i) = SE3_AdT(invT(:, :, i))*V(:, i-1) + dq(i)*L(:, i);
    dV(:, i) = SE3_AdT(invT(:, :, i))*dV(:, i-1) + dq(i)*SE3_adxi(V(:, i))*L(:, i) + ddq(i)*L(:, i);
end

% Backward iterations
F(:, n) = -SE3_AdT(T0(:, :, n+1))*fext + I(:, :, n)*dV(:, n) - SE3_adxi(V(:, n)).'*I(:, :, n)*V(:, n);
tau(n) = L(:, n).'*F(:, n);
for i = n-1:-1:1
    F(:, i) = SE3_AdT(invT(:, :, i+1)).'*F(:, i+1) + I(:, :, i)*dV(:, i) - SE3_adxi(V(:, i)).'*I(:, :, i)*V(:, i);
    tau(i) = L(:, i).'*F(:, i);
end
end