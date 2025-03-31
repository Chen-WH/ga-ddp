function [p2tau_pqpq, p2tau_pqpdq, p2tau_pdqpdq, p2tau_pqpddq] = GA_idyn_so(rbt, q, dq, ddq, fext)
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
M(:, 1) = GA_prodM(M0(:, 1), GA_motor(L(:, 1), q(1), type(1)));
Mt(:, 1) = GA_Rev(M(:, 1));
V(:, 1) = dq(1)*L(:, 1);
dV(:, 1) = PGA_Frame2([-rbt.g; 0; 0; 0], Mt(:, 1)) + ddq(1)*L(:, 1);

pV_pdq(:, 1, 1) = L(:, 1);
pdV_pddq(:, 1, 1) = L(:, 1);

for k = 2:n
    % Kinect forward pass
    M(:, k) = GA_prodM(M0(:, k), GA_motor(L(:, k), q(k), type(k)));
    Mt(:, k) = GA_Rev(M(:, k));
    AdM = GA_AdM(Mt(:, k));
    V(:, k) = PGA_Frame2(V(:, k-1), Mt(:, k)) + dq(k)*L(:, k);
    dV(:, k) = PGA_Frame2(dV(:, k-1), Mt(:, k)) - dq(k)*PGA_Lie2(L(:, k), V(:, k)) + ddq(k)*L(:, k);
    adL = GA_adB(L(:, k));

    % first-order pass
    pV_pq(:, 1:k, k) = AdM*pV_pq(:, 1:k, k-1);
    pV_pq(:, k, k) = pV_pq(:, k, k) - PGA_Lie2(L(:, k), V(:, k)); % i == k case

    pV_pdq(:, 1:k, k) = AdM*pV_pdq(:, 1:k, k-1);
    pV_pdq(:, k, k) = pV_pdq(:, k, k) + L(:, k); % i == k case

    pdV_pq(:, 1:k, k) = AdM*pdV_pq(:, 1:k, k-1) - dq(k)*adL*pV_pq(:, 1:k, k);
    pdV_pq(:, k, k) = pdV_pq(:, k, k) - PGA_Lie2(L(:, k), PGA_Frame2(dV(:, k-1), Mt(:, k))); % i == k case

    pdV_pdq(:, 1:k, k) = AdM*pdV_pdq(:, 1:k, k-1) - dq(k)*adL*pV_pdq(:, 1:k, k);
    pdV_pdq(:, k, k) = pdV_pdq(:, k, k) - PGA_Lie2(L(:, k), V(:, k)); % i == k case
    
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
    p2dV_pqpq(:, k, k, k)  = p2dV_pqpq(:, k, k, k) + PGA_Lie2(L(:, k), PGA_Lie2(L(:, k), PGA_Frame2(dV(:, k-1), Mt(:, k)))); % i == j == k case

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
F(:, n) = PGA_Frame2(fext, M0(:, n+1)) + I(:, :, n)*dV(:, n) - PGA_Lie2(I(:, :, n)*V(:, n), V(:, n));
adIV = GA_adB(I(:, :, n)*V(:, n));
adV = GA_adB(V(:, n));

pF_pq(:, :, n) = I(:, :, n)*pdV_pq(:, :, n) - adIV*pV_pq(:, :, n) + adV*I(:, :, n)*pV_pq(:, :, n);
pF_pdq(:, :, n) = I(:, :, n)*pdV_pdq(:, :, n) - adIV*pV_pdq(:, :, n) + adV*I(:, :, n)*pV_pdq(:, :, n);
pF_pddq(:, :, n) = I(:, :, n)*pdV_pddq(:, :, n);

tmp = pagemtimes(GA_adB(I(:, :, n)*pV_pq(:, :, n)), pV_pq(:, :, n));
p2F_pqpq(:, :, :, n) = reshape(I(:, :, n)*reshape(p2dV_pqpq(:, :, :, n), 6, n^2) + adV*I(:, :, n)*reshape(p2V_pqpq(:, :, :, n), 6, n^2) - adIV*reshape(p2V_pqpq(:, :, :, n), 6, n^2), 6, n, n) - tmp - permute(tmp, [1, 3, 2]);
p2F_pqpdq(:, :, :, n) = reshape(I(:, :, n)*reshape(p2dV_pqpdq(:, :, :, n), 6, n^2) + adV*I(:, :, n)*reshape(p2V_pqpdq(:, :, :, n), 6, n^2) - adIV*reshape(p2V_pqpdq(:, :, :, n), 6, n^2), 6, n, n)...
    - pagemtimes(GA_adB(I(:, :, n)*pV_pdq(:, :, n)), pV_pq(:, :, n)) - permute(pagemtimes(GA_adB(I(:, :, n)*pV_pq(:, :, n)), pV_pdq(:, :, n)), [1, 3, 2]);
tmp = pagemtimes(GA_adB(I(:, :, n)*pV_pdq(:, :, n)), pV_pdq(:, :, n));
p2F_pdqpdq(:, :, :, n) = reshape(I(:, :, n)*reshape(p2dV_pdqpdq(:, :, :, n), 6, n^2), 6, n, n) - tmp - permute(tmp, [1, 3, 2]);
p2F_pqpddq(:, :, :, n) = reshape(I(:, :, n)*reshape(p2dV_pqpddq(:, :, :, n), 6, n^2), 6, n, n);

p2tau_pqpq(n, :, :) = reshape(GA_metric(L(:, n))*reshape(p2F_pqpq(:, :, :, n), 6, n^2), 1, n, n);
p2tau_pqpdq(n, :, :) = reshape(GA_metric(L(:, n))*reshape(p2F_pqpdq(:, :, :, n), 6, n^2), 1, n, n);
p2tau_pdqpdq(n, :, :) = reshape(GA_metric(L(:, n))*reshape(p2F_pdqpdq(:, :, :, n), 6, n^2), 1, n, n);
p2tau_pqpddq(n, :, :) = reshape(GA_metric(L(:, n))*reshape(p2F_pqpddq(:, :, :, n), 6, n^2), 1, n, n);

for k = n-1:-1:1
    % ID pass
    AdM = GA_AdM(M(:, k+1));
    F(:, k) = PGA_Frame2(F(:, k+1), M(:, k+1)) + I(:, :, k)*dV(:, k) - PGA_Lie2(I(:, :, k)*V(:, k), V(:, k));
    adL = GA_adB(L(:, k+1));
    adIV = GA_adB(I(:, :, k)*V(:, k));
    adV = GA_adB(V(:, k));

    % first-order pass
    pF_pq(:, :, k) = AdM*pF_pq(:, :, k+1) + I(:, :, k)*pdV_pq(:, :, k) - adIV*pV_pq(:, :, k) + adV*I(:, :, k)*pV_pq(:, :, k);
    pF_pq(:, k+1, k) = pF_pq(:, k+1, k) + PGA_Frame2(PGA_Lie2(L(:, k+1), F(:, k+1)), M(:, k+1)); % i == k+1 case
    pF_pdq(:, :, k) = AdM*pF_pdq(:, :, k+1) + I(:, :, k)*pdV_pdq(:, :, k) - adIV*pV_pdq(:, :, k) + adV*I(:, :, k)*pV_pdq(:, :, k);
    pF_pddq(:, :, k) = AdM*pF_pddq(:, :, k+1) + I(:, :, k)*pdV_pddq(:, :, k);

    % second-order pass
    IpV_pq = GA_adB(I(:, :, k)*pV_pq(:, :, k));
    IpV_pdq = GA_adB(I(:, :, k)*pV_pdq(:, :, k));

    tmp = pagemtimes(IpV_pq, pV_pq(:, :, k));
    p2F_pqpq(:, :, :, k) = reshape(AdM*reshape(p2F_pqpq(:, :, :, k+1), 6, n^2) + I(:, :, k)*reshape(p2dV_pqpq(:, :, :, k), 6, n^2) + adV*I(:, :, k)*reshape(p2V_pqpq(:, :, :, k), 6, n^2) - adIV*reshape(p2V_pqpq(:, :, :, k), 6, n^2), 6, n, n) - tmp - permute(tmp, [1, 3, 2]);
    p2F_pqpq(:, k+1, :, k) = p2F_pqpq(:, k+1, :, k) + reshape(AdM*adL*pF_pq(:, :, k+1), 6, 1, n); % i == k+1 case
    p2F_pqpq(:, :, k+1, k) = p2F_pqpq(:, :, k+1, k) + AdM*adL*pF_pq(:, :, k+1); % j == k+1 case
    p2F_pqpq(:, k+1, k+1, k) = p2F_pqpq(:, k+1, k+1, k) + PGA_Frame2(PGA_Lie2(L(:, k+1), PGA_Lie2(L(:, k+1), F(:, k+1))), M(:, k+1)); % i == j == k+1 case
    
    p2F_pqpdq(:, :, :, k) = reshape(AdM*reshape(p2F_pqpdq(:, :, :, k+1), 6, n^2) + I(:, :, k)*reshape(p2dV_pqpdq(:, :, :, k), 6, n^2) + adV*I(:, :, k)*reshape(p2V_pqpdq(:, :, :, k), 6, n^2) - adIV*reshape(p2V_pqpdq(:, :, :, k), 6, n^2), 6, n, n)...
        - pagemtimes(IpV_pdq, pV_pq(:, :, k)) - permute(pagemtimes(IpV_pq, pV_pdq(:, :, k)), [1, 3, 2]);
    p2F_pqpdq(:, k+1, :, k) = p2F_pqpdq(:, k+1, :, k) + reshape(AdM*adL*pF_pdq(:, :, k+1), 6, 1, n); % i == k+1 case

    tmp = pagemtimes(IpV_pdq, pV_pdq(:, :, k));
    p2F_pdqpdq(:, :, :, k) = reshape(AdM*reshape(p2F_pdqpdq(:, :, :, k+1), 6, n^2) + I(:, :, k)*reshape(p2dV_pdqpdq(:, :, :, k), 6, n^2), 6, n, n) - tmp - permute(tmp, [1, 3, 2]);
    
    p2F_pqpddq(:, :, :, k) = reshape(AdM*reshape(p2F_pqpddq(:, :, :, k+1), 6, n^2) + I(:, :, k)*reshape(p2dV_pqpddq(:, :, :, k), 6, n^2), 6, n, n);
    p2F_pqpddq(:, k+1, :, k) = p2F_pqpddq(:, k+1, :, k) + reshape(AdM*adL*pF_pddq(:, :, k+1), 6, 1, n); % i == k+1 case
    
    p2tau_pqpq(k, :, :) = reshape(GA_metric(L(:, k))*reshape(p2F_pqpq(:, :, :, k), 6, n^2), 1, n, n);
    p2tau_pqpdq(k, :, :) = reshape(GA_metric(L(:, k))*reshape(p2F_pqpdq(:, :, :, k), 6, n^2), 1, n, n);
    p2tau_pdqpdq(k, :, :) = reshape(GA_metric(L(:, k))*reshape(p2F_pdqpdq(:, :, :, k), 6, n^2), 1, n, n);
    p2tau_pqpddq(k, :, :) = reshape(GA_metric(L(:, k))*reshape(p2F_pqpddq(:, :, :, k), 6, n^2), 1, n, n);
end

end