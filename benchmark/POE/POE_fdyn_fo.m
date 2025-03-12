function [pddq_pq, pddq_pdq, pddq_ptau] = POE_fdyn_fo(rbt, q, dq, tau, fext)
if nargin == 4
    fext = zeros(6, 1);
end
dual = @(x) [x(4:6, :); x(1:3, :)];
ddq = POE_fdyn(rbt, q, dq, tau);

n = rbt.n;
T0 = rbt.Tj;
L = rbt.Lj;
I = rbt.Ij;
type = rbt.type;

% init
invT = zeros(4, 4, n);
V = zeros(6, n);
pV_pq = zeros(6, n, n);
pV_pdq = zeros(6, n, n);
dV = zeros(6, n);
pdV_pq = zeros(6, n, n);
pdV_pdq = zeros(6, n, n);
pdV_pddq = zeros(6, n, n);
F = zeros(6, n);
pF_pq = zeros(6, n, n);
pF_pdq = zeros(6, n, n);
pF_pddq = zeros(6, n, n);
ptau_pq = zeros(n, n);
ptau_pdq = zeros(n, n);
ptau_pddq = zeros(n, n);

% Forward iterations
invT(:, :, 1) = SE3_inv(T0(:, :, 1)*POE_motor(L(:, 1), q(1), type(1)));
V(:, 1) = dq(1)*L(:, 1);
dV(:, 1) = SE3_AdT(invT(:, :, 1))*[0; 0; 0; -rbt.g] + ddq(1)*L(:, 1);

pV_pdq(:, 1, 1)= L(:, 1);
pdV_pddq(:, 1, 1) = L(:, 1);

for k = 2:n
    % Kinect pass
    invT(:, :, k) = SE3_inv(T0(:, :, k)*POE_motor(L(:, k), q(k), type(k)));
    AdM = SE3_AdT(invT(:, :, k));
    adL = SE3_adxi(L(:, k));
    V(:, k) = AdM*V(:, k-1) + dq(k)*L(:, k);
    dV(:, k) = AdM*dV(:, k-1) - dq(k)*adL*V(:, k) + ddq(k)*L(:, k);

    % first-order pass
    pV_pq(:, 1:k, k) = AdM*pV_pq(:, 1:k, k-1);
    pV_pq(:, k, k) = pV_pq(:, k, k) - adL*V(:, k); % i == k case

    pV_pdq(:, 1:k, k) = AdM*pV_pdq(:, 1:k, k-1);
    pV_pdq(:, k, k) = pV_pdq(:, k, k) + L(:, k); % i == k case

    pdV_pq(:, 1:k, k) = AdM*pdV_pq(:, 1:k, k-1) - dq(k)*adL*pV_pq(:, 1:k, k);
    pdV_pq(:, k, k) = pdV_pq(:, k, k) - adL*AdM*dV(:, k-1); % i == k case

    pdV_pdq(:, 1:k, k) = AdM*pdV_pdq(:, 1:k, k-1) - dq(k)*adL*pV_pdq(:, 1:k, k);
    pdV_pdq(:, k, k) = pdV_pdq(:, k, k) - adL*V(:, k); % i == k case

    pdV_pddq(:, 1:k, k) = AdM*pdV_pddq(:, 1:k, k-1);
    pdV_pddq(:, k, k) = pdV_pddq(:, k, k) + L(:, k); % i == k case
end

% Backward iterations
F(:, n) = -SE3_AdT(T0(:, :, n+1))*fext + I(:, :, n)*dV(:, n) - SE3_adxi(V(:, n)).'*I(:, :, n)*V(:, n);
adIV = SE3_adxi(dual(I(:, :, n)*V(:, n)));
adV = SE3_adxi(V(:, n)).';

pF_pq(:, :, n) = I(:, :, n)*pdV_pq(:, :, n) - dual(adIV*pV_pq(:, :, n)) - adV*I(:, :, n)*pV_pq(:, :, n);
pF_pdq(:, :, n) = I(:, :, n)*pdV_pdq(:, :, n) - dual(adIV*pV_pdq(:, :, n)) - adV*I(:, :, n)*pV_pdq(:, :, n);
pF_pddq(:, :, n) = I(:, :, n)*pdV_pddq(:, :, n);
ptau_pq(n, :) = L(:, k).'*pF_pq(:, :, n);
ptau_pdq(n, :) = L(:, k).'*pF_pdq(:, :, n);
ptau_pddq(n, :) = L(:, k).'*pF_pddq(:, :, n);

for k = n-1:-1:1
    % ID pass
    AdM = SE3_AdT(invT(:, :, k+1)).';
    adIV = SE3_adxi(dual(I(:, :, k)*V(:, k)));
    adV = SE3_adxi(V(:, k)).';
    F(:, k) = AdM*F(:, k+1) + I(:, :, k)*dV(:, k) - adV*I(:, :, k)*V(:, k);

    % first-order pass
    pF_pq(:, :, k) = AdM*pF_pq(:, :, k+1) + I(:, :, k)*pdV_pq(:, :, k) - dual(adIV*pV_pq(:, :, k)) - adV*I(:, :, k)*pV_pq(:, :, k);
    pF_pq(:, k+1, k) = pF_pq(:, k+1, k) + AdM*SE3_adxi(-L(:, k+1)).'*F(:, k+1); % i == k+1 case
    pF_pdq(:, :, k) = AdM*pF_pdq(:, :, k+1) + I(:, :, k)*pdV_pdq(:, :, k) - dual(adIV*pV_pdq(:, :, k)) - adV*I(:, :, k)*pV_pdq(:, :, k);
    pF_pddq(:, :, k) = AdM*pF_pddq(:, :, k+1) + I(:, :, k)*pdV_pddq(:, :, k);

    ptau_pq(k, :) = L(:, k).'*pF_pq(:, :, k);
    ptau_pdq(k, :) = L(:, k).'*pF_pdq(:, :, k);
    ptau_pddq(k, :) = L(:, k).'*pF_pddq(:, :, k);
end

pddq_ptau = inv(ptau_pddq);
pddq_pq = -pddq_ptau*ptau_pq;
pddq_pdq = -pddq_ptau*ptau_pdq;
end