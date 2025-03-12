function [ptau_pq, ptau_pdq, ptau_pddq] = GA_idyn_fo(rbt, q, dq, ddq, fext)
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
end

% Backward iterations
F(:, n) = PGA_Frame2(fext, M0(:, n+1)) + I(:, :, n)*dV(:, n) - PGA_Lie2(I(:, :, n)*V(:, n), V(:, n));
adIV = GA_adB(I(:, :, n)*V(:, n));
adV = GA_adB(V(:, n));

pF_pq(:, :, n) = I(:, :, n)*pdV_pq(:, :, n) - adIV*pV_pq(:, :, n) + adV*I(:, :, n)*pV_pq(:, :, n);
pF_pdq(:, :, n) = I(:, :, n)*pdV_pdq(:, :, n) - adIV*pV_pdq(:, :, n) + adV*I(:, :, n)*pV_pdq(:, :, n);
pF_pddq(:, :, n) = I(:, :, n)*pdV_pddq(:, :, n);
ptau_pq(n, :) = GA_metric(L(:, n))*pF_pq(:, :, n);
ptau_pdq(n, :) = GA_metric(L(:, n))*pF_pdq(:, :, n);
ptau_pddq(n, :) = GA_metric(L(:, n))*pF_pddq(:, :, n);

for k = n-1:-1:1
    % Kinect backward pass
    AdM = GA_AdM(M(:, k+1));
    F(:, k) = PGA_Frame2(F(:, k+1), M(:, k+1)) + I(:, :, k)*dV(:, k) - PGA_Lie2(I(:, :, k)*V(:, k), V(:, k));
    adIV = GA_adB(I(:, :, k)*V(:, k));
    adV = GA_adB(V(:, k));

    % first-order pass
    pF_pq(:, :, k) = AdM*pF_pq(:, :, k+1) + I(:, :, k)*pdV_pq(:, :, k) - adIV*pV_pq(:, :, k) + adV*I(:, :, k)*pV_pq(:, :, k);
    pF_pq(:, k+1, k) = pF_pq(:, k+1, k) + PGA_Frame2(PGA_Lie2(L(:, k+1), F(:, k+1)), M(:, k+1)); % i == k+1 case
    pF_pdq(:, :, k) = AdM*pF_pdq(:, :, k+1) + I(:, :, k)*pdV_pdq(:, :, k) - adIV*pV_pdq(:, :, k) + adV*I(:, :, k)*pV_pdq(:, :, k);
    pF_pddq(:, :, k) = AdM*pF_pddq(:, :, k+1) + I(:, :, k)*pdV_pddq(:, :, k);

    ptau_pq(k, :) = GA_metric(L(:, k))*pF_pq(:, :, k);
    ptau_pdq(k, :) = GA_metric(L(:, k))*pF_pdq(:, :, k);
    ptau_pddq(k, :) = GA_metric(L(:, k))*pF_pddq(:, :, k);
end
end