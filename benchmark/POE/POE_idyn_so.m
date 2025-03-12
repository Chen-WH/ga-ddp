function [p2tau_pqpq, p2tau_pqpdq, p2tau_pdqpdq, p2tau_pqpddq] = POE_idyn_so(rbt, q, dq, ddq, fext)
if nargin == 4
    fext = zeros(6, 1);
end
dual = @(x) [x(4:6, :, :); x(1:3, :, :)];
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
p2V_pqpq = zeros(6, n, n, n);
p2V_pqpdq = zeros(6, n, n, n);
dV = zeros(6, n);
pdV_pq = zeros(6, n, n);
pdV_pdq = zeros(6, n, n);
pdV_pddq = zeros(6, n, n);
p2dV_pqpq = zeros(6, n, n, n);
p2dV_pqpdq = zeros(6, n, n, n);
p2dV_pdqpdq = zeros(6, n, n, n);
p2dV_pqpddq = zeros(6, n, n, n);
F = zeros(6, n);
pF_pq = zeros(6, n, n);
pF_pdq = zeros(6, n, n);
pF_pddq = zeros(6, n, n);
p2F_pqpq = zeros(6, n, n, n);
p2F_pqpdq = zeros(6, n, n, n);
p2F_pdqpdq = zeros(6, n, n, n);
p2F_pqpddq = zeros(6, n, n, n);
p2tau_pqpq = zeros(n, n, n);
p2tau_pqpdq = zeros(n, n, n);
p2tau_pdqpdq = zeros(n, n, n);
p2tau_pqpddq = zeros(n, n, n);

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
    V(:, k) = SE3_AdT(invT(:, :, k))*V(:, k-1) + dq(k)*L(:, k);
    dV(:, k) = SE3_AdT(invT(:, :, k))*dV(:, k-1) - dq(k)*adL*V(:, k) + ddq(k)*L(:, k);

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

    % second-order pass
    p2V_pqpq(:, 1:k, 1:k, k) = reshape(AdM*reshape(p2V_pqpq(:, 1:k, 1:k, k-1), 6, k^2), 6, k, k);
    p2V_pqpq(:, k, 1:k, k)  = p2V_pqpq(:, k, 1:k, k) - reshape(adL*pV_pq(:, 1:k, k), 6, 1, k); % i == k case
    p2V_pqpq(:, 1:k, k, k)  = p2V_pqpq(:, 1:k, k, k) - adL*AdM*pV_pq(:, 1:k, k-1); % j == k case

    p2V_pqpdq(:, 1:k, 1:k, k) = reshape(AdM*reshape(p2V_pqpdq(:, 1:k, 1:k, k-1), 6, k^2), 6, k, k);
    p2V_pqpdq(:, k, 1:k, k)  = p2V_pqpdq(:, k, 1:k, k) - reshape(adL*pV_pdq(:, 1:k, k), 6, 1, k); % i == k case

    p2dV_pqpq(:, 1:k, 1:k, k) = reshape(AdM*reshape(p2dV_pqpq(:, 1:k, 1:k, k-1), 6, k^2), 6, k, k) - reshape(dq(k)*adL*reshape(p2V_pqpq(:, 1:k, 1:k, k), 6, k^2), 6, k, k);
    p2dV_pqpq(:, k, 1:k, k)  = p2dV_pqpq(:, k, 1:k, k) - reshape(adL*AdM*pdV_pq(:, 1:k, k-1), 6, 1, k); % i == k case
    p2dV_pqpq(:, 1:k, k, k)  = p2dV_pqpq(:, 1:k, k, k) - adL*AdM*pdV_pq(:, 1:k, k-1); % j == k case
    p2dV_pqpq(:, k, k, k)  = p2dV_pqpq(:, k, k, k) + adL*adL*AdM*dV(:, k-1); % i == j == k case

    p2dV_pqpdq(:, 1:k, 1:k, k) = reshape(AdM*reshape(p2dV_pqpdq(:, 1:k, 1:k, k-1), 6, k^2), 6, k, k) - reshape(dq(k)*adL*reshape(p2V_pqpdq(:, 1:k, 1:k, k), 6, k^2), 6, k, k);
    p2dV_pqpdq(:, k, 1:k, k)  = p2dV_pqpdq(:, k, 1:k, k) - reshape(adL*AdM*pdV_pdq(:, 1:k, k-1), 6, 1, k); % i == k case
    p2dV_pqpdq(:, 1:k, k, k)  = p2dV_pqpdq(:, 1:k, k, k) - adL*pV_pq(:, 1:k, k); % j == k case

    p2dV_pdqpdq(:, 1:k, 1:k, k) = reshape(AdM*reshape(p2dV_pdqpdq(:, 1:k, 1:k, k-1), 6, k^2), 6, k, k);
    p2dV_pdqpdq(:, k, 1:k, k)  = p2dV_pdqpdq(:, k, 1:k, k) - reshape(adL*pV_pdq(:, 1:k, k), 6, 1, k); % i == k case
    p2dV_pdqpdq(:, 1:k, k, k)  = p2dV_pdqpdq(:, 1:k, k, k) - adL*pV_pdq(:, 1:k, k); % j == k case

    p2dV_pqpddq(:, 1:k, 1:k, k) = reshape(AdM*reshape(p2dV_pqpddq(:, 1:k, 1:k, k-1), 6, k^2), 6, k, k);
    p2dV_pqpddq(:, k, 1:k, k)  = p2dV_pqpddq(:, k, 1:k, k) - reshape(adL*AdM*pdV_pddq(:, 1:k, k-1), 6, 1, k); % i == k case
end

% Backward iterations
F(:, n) = -SE3_AdT(T0(:, :, n+1))*fext + I(:, :, n)*dV(:, n) - SE3_adxi(V(:, n)).'*I(:, :, n)*V(:, n);
adIV = SE3_adxi(dual(I(:, :, n)*V(:, n)));
adV = SE3_adxi(V(:, n)).';

pF_pq(:, :, n) = I(:, :, n)*pdV_pq(:, :, n) - dual(adIV*pV_pq(:, :, n)) - adV*I(:, :, n)*pV_pq(:, :, n);
pF_pdq(:, :, n) = I(:, :, n)*pdV_pdq(:, :, n) - dual(adIV*pV_pdq(:, :, n)) - adV*I(:, :, n)*pV_pdq(:, :, n);
pF_pddq(:, :, n) = I(:, :, n)*pdV_pddq(:, :, n);

tmp = dual(pagemtimes(SE3_adxi(dual(I(:, :, n)*pV_pq(:, :, n))), pV_pq(:, :, n)));
p2F_pqpq(:, :, :, n) = reshape(I(:, :, n)*reshape(p2dV_pqpq(:, :, :, n), 6, n^2) - adV*I(:, :, n)*reshape(p2V_pqpq(:, :, :, n), 6, n^2) - dual(adIV*reshape(p2V_pqpq(:, :, :, n), 6, n^2)), 6, n, n) - tmp - permute(tmp, [1, 3, 2]);
p2F_pqpdq(:, :, :, n) = reshape(I(:, :, n)*reshape(p2dV_pqpdq(:, :, :, n), 6, n^2) - adV*I(:, :, n)*reshape(p2V_pqpdq(:, :, :, n), 6, n^2) - dual(adIV*reshape(p2V_pqpdq(:, :, :, n), 6, n^2)), 6, n, n)...
    - dual(pagemtimes(SE3_adxi(dual(I(:, :, n)*pV_pdq(:, :, n))), pV_pq(:, :, n))) - permute(dual(pagemtimes(SE3_adxi(dual(I(:, :, n)*pV_pq(:, :, n))), pV_pdq(:, :, n))), [1, 3, 2]);
tmp = dual(pagemtimes(SE3_adxi(dual(I(:, :, n)*pV_pdq(:, :, n))), pV_pdq(:, :, n)));
p2F_pdqpdq(:, :, :, n) = reshape(I(:, :, n)*reshape(p2dV_pdqpdq(:, :, :, n), 6, n^2), 6, n, n) - tmp - permute(tmp, [1, 3, 2]);
p2F_pqpddq(:, :, :, n) = reshape(I(:, :, n)*reshape(p2dV_pqpddq(:, :, :, n), 6, n^2), 6, n, n);

p2tau_pqpq(n, :, :) = reshape(L(:, k).'*reshape(p2F_pqpq(:, :, :, n), 6, n^2), 1, n, n);
p2tau_pqpdq(n, :, :) = reshape(L(:, k).'*reshape(p2F_pqpdq(:, :, :, n), 6, n^2), 1, n, n);
p2tau_pdqpdq(n, :, :) = reshape(L(:, k).'*reshape(p2F_pdqpdq(:, :, :, n), 6, n^2), 1, n, n);
p2tau_pqpddq(n, :, :) = reshape(L(:, k).'*reshape(p2F_pqpddq(:, :, :, n), 6, n^2), 1, n, n);

for k = n-1:-1:1
    % ID pass
    AdM = SE3_AdT(invT(:, :, k+1)).';
    adIV = SE3_adxi(dual(I(:, :, k)*V(:, k)));
    adV = SE3_adxi(V(:, k)).';
    F(:, k) = AdM*F(:, k+1) + I(:, :, k)*dV(:, k) - adV*I(:, :, k)*V(:, k);

    % first-order pass
    pF_pq(:, :, k) = AdM*pF_pq(:, :, k+1) + I(:, :, k)*pdV_pq(:, :, k) - dual(adIV*pV_pq(:, :, k)) - adV*I(:, :, k)*pV_pq(:, :, k);
    pF_pq(:, k+1, k) = pF_pq(:, k+1, k) + AdM*adL*F(:, k+1); % i == k+1 case
    pF_pdq(:, :, k) = AdM*pF_pdq(:, :, k+1) + I(:, :, k)*pdV_pdq(:, :, k) - dual(adIV*pV_pdq(:, :, k)) - adV*I(:, :, k)*pV_pdq(:, :, k);
    pF_pddq(:, :, k) = AdM*pF_pddq(:, :, k+1) + I(:, :, k)*pdV_pddq(:, :, k);

    % second-order pass
    IpV_pq = SE3_adxi(dual(I(:, :, k)*pV_pq(:, :, k)));
    IpV_pdq = SE3_adxi(dual(I(:, :, k)*pV_pdq(:, :, k)));

    tmp = dual(pagemtimes(IpV_pq, pV_pq(:, :, k)));
    p2F_pqpq(:, :, :, k) = reshape(AdM*reshape(p2F_pqpq(:, :, :, k+1), 6, n^2) + I(:, :, k)*reshape(p2dV_pqpq(:, :, :, k), 6, n^2) - adV*I(:, :, k)*reshape(p2V_pqpq(:, :, :, k), 6, n^2) - dual(adIV*reshape(p2V_pqpq(:, :, :, k), 6, n^2)), 6, n, n) - tmp - permute(tmp, [1, 3, 2]);
    p2F_pqpq(:, k+1, :, k) = p2F_pqpq(:, k+1, :, k) + reshape(AdM*adL*pF_pq(:, :, k+1), 6, 1, n); % i == k+1 case
    p2F_pqpq(:, :, k+1, k) = p2F_pqpq(:, :, k+1, k) + AdM*adL*pF_pq(:, :, k+1); % j == k+1 case
    p2F_pqpq(:, k+1, k+1, k) = p2F_pqpq(:, k+1, k+1, k) + AdM*adL*adL*F(:, k+1); % i == j == k+1 case
    
    p2F_pqpdq(:, :, :, k) = reshape(AdM*reshape(p2F_pqpdq(:, :, :, k+1), 6, n^2) + I(:, :, k)*reshape(p2dV_pqpdq(:, :, :, k), 6, n^2) - adV*I(:, :, k)*reshape(p2V_pqpdq(:, :, :, k), 6, n^2) - dual(adIV*reshape(p2V_pqpdq(:, :, :, k), 6, n^2)), 6, n, n)...
        - dual(pagemtimes(IpV_pdq, pV_pq(:, :, k))) - permute(dual(pagemtimes(IpV_pq, pV_pdq(:, :, k))), [1, 3, 2]);
    p2F_pqpdq(:, k+1, :, k) = p2F_pqpdq(:, k+1, :, k) + reshape(AdM*adL*pF_pdq(:, :, k+1), 6, 1, n); % i == k+1 case

    tmp = dual(pagemtimes(IpV_pdq, pV_pdq(:, :, k)));
    p2F_pdqpdq(:, :, :, k) = reshape(AdM*reshape(p2F_pdqpdq(:, :, :, k+1), 6, n^2) + I(:, :, k)*reshape(p2dV_pdqpdq(:, :, :, k), 6, n^2), 6, n, n) - tmp - permute(tmp, [1, 3, 2]);
    
    p2F_pqpddq(:, :, :, k) = reshape(AdM*reshape(p2F_pqpddq(:, :, :, k+1), 6, n^2) + I(:, :, k)*reshape(p2dV_pqpddq(:, :, :, k), 6, n^2), 6, n, n);
    p2F_pqpddq(:, k+1, :, k) = p2F_pqpddq(:, k+1, :, k) + reshape(AdM*adL*pF_pddq(:, :, k+1), 6, 1, n); % i == k+1 case
    
    p2tau_pqpq(k, :, :) = reshape(L(:, k).'*reshape(p2F_pqpq(:, :, :, k), 6, n^2), 1, n, n);
    p2tau_pqpdq(k, :, :) = reshape(L(:, k).'*reshape(p2F_pqpdq(:, :, :, k), 6, n^2), 1, n, n);
    p2tau_pdqpdq(k, :, :) = reshape(L(:, k).'*reshape(p2F_pdqpdq(:, :, :, k), 6, n^2), 1, n, n);
    p2tau_pqpddq(k, :, :) = reshape(L(:, k).'*reshape(p2F_pqpddq(:, :, :, k), 6, n^2), 1, n, n);
end
end