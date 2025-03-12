clc; clear; close all;
rbt_ga = GA_panda();
rbt_poe = POE_panda();
rbt_rst = importrobot(rbt_ga.urdf, DataFormat='column');
rbt_rst.Gravity = rbt_ga.g;
q = rand(rbt_ga.n, 1)*2*pi - pi;
dq = rand(rbt_ga.n, 1)*10 - 5;
ddq = rand(rbt_ga.n, 1)*10 - 5;
% q = [-1.2006; 1.4207; 1.7773; 1.2176; -3.0800; 2.1565];
% dq = [4.2233; 2.7095; -4.5734; -1.2181; 2.0434; 2.2951];
% ddq = [-2.7572; -2.3095; 1.7303; -0.2251; 1.2372; -2.6356];
tau = inverseDynamics(rbt_rst, q, dq, ddq);

%% GA
[p2tau_pqpq_ga, p2tau_pqpdq_ga, p2tau_pdqpdq_ga, p2tau_pqpddq_ga] = GA_idyn_so(rbt_ga, q, dq, ddq);
[p2ddq_pqpq_ga, p2ddq_pdqpq_ga, p2ddq_pdqpdq_ga, p2ddq_ptaupq_ga] = GA_fdyn_so(rbt_ga, q, dq, tau);

%% POE
[p2tau_pqpq_poe, p2tau_pqpdq_poe, p2tau_pdqpdq_poe, p2tau_pqpddq_poe] = POE_idyn_so(rbt_poe, q, dq, ddq);
[p2ddq_pqpq_poe, p2ddq_pdqpq_poe, p2ddq_pdqpdq_poe, p2ddq_ptaupq_poe] = POE_fdyn_so(rbt_poe, q, dq, tau);

%% Finite difference
p2tau_pqpq_rst = zeros(rbt_ga.n, rbt_ga.n, rbt_ga.n);
p2tau_pqpdq_rst = zeros(rbt_ga.n, rbt_ga.n, rbt_ga.n);
p2tau_pdqpdq_rst = zeros(rbt_ga.n, rbt_ga.n, rbt_ga.n);
p2tau_pqpddq_rst = zeros(rbt_ga.n, rbt_ga.n, rbt_ga.n);
p2ddq_pqpq_rst = zeros(rbt_ga.n, rbt_ga.n, rbt_ga.n);
p2ddq_pdqpq_rst = zeros(rbt_ga.n, rbt_ga.n, rbt_ga.n);
p2ddq_pdqpdq_rst = zeros(rbt_ga.n, rbt_ga.n, rbt_ga.n);
p2ddq_ptaupq_rst = zeros(rbt_ga.n, rbt_ga.n, rbt_ga.n);
step = 1e-7;
for j = 1:rbt_ga.n
    [ptau_pq, ptau_pdq, ~] = GA_idyn_fo(rbt_ga, q, dq, ddq);
    [pddq_pq, pddq_pdq, pddq_ptau] = GA_fdyn_fo(rbt_ga, q, dq, tau);
    q_tmp = q;
    dq_tmp = dq;
    ddq_tmp = ddq;
    tau_tmp = tau;
    q_tmp(j) = q_tmp(j) + step;
    dq_tmp(j) = dq_tmp(j) + step;
    ddq_tmp(j) = ddq_tmp(j) + step;
    % ID-SO
    [ptau_pq_tmp, ~, ~] = GA_idyn_fo(rbt_ga, q_tmp, dq, ddq);
    p2tau_pqpq_rst(:, :, j) = (ptau_pq_tmp - ptau_pq)/step;
    [ptau_pq_tmp, ptau_pdq_tmp, ~] = GA_idyn_fo(rbt_ga, q, dq_tmp, ddq);
    p2tau_pqpdq_rst(:, :, j) = (ptau_pq_tmp - ptau_pq)/step;
    p2tau_pdqpdq_rst(:, :, j) = (ptau_pdq_tmp - ptau_pdq)/step;
    [ptau_pq_tmp, ~, ~] = GA_idyn_fo(rbt_ga, q, dq, ddq_tmp);
    p2tau_pqpddq_rst(:, :, j) = (ptau_pq_tmp - ptau_pq)/step;
    % FD-SO
    [pddq_pq_tmp, pddq_pdq_tmp, pddq_ptau_tmp] = GA_fdyn_fo(rbt_ga, q_tmp, dq, tau);
    p2ddq_pqpq_rst(:, :, j) = (pddq_pq_tmp - pddq_pq)/step;
    p2ddq_pdqpq_rst(:, :, j) = (pddq_pdq_tmp - pddq_pdq)/step;
    p2ddq_ptaupq_rst(:, :, j) = (pddq_ptau_tmp - pddq_ptau)/step;
    [~, pddq_pdq_tmp, ~] = GA_fdyn_fo(rbt_ga, q, dq_tmp, tau);
    p2ddq_pdqpdq_rst(:, :, j) = (pddq_pdq_tmp - pddq_pdq)/step;
end

%% Result Validation
norm3d = @(A) sqrt(sum(A(:).^2));
norm3d(p2tau_pqpq_ga - p2tau_pqpq_rst)
norm3d(p2tau_pqpdq_ga - p2tau_pqpdq_rst)
norm3d(p2tau_pdqpdq_ga - p2tau_pdqpdq_rst)
norm3d(p2tau_pqpddq_ga - p2tau_pqpddq_rst)

norm3d(p2tau_pqpq_poe - p2tau_pqpq_rst)
norm3d(p2tau_pqpdq_poe - p2tau_pqpdq_rst)
norm3d(p2tau_pdqpdq_poe - p2tau_pdqpdq_rst)
norm3d(p2tau_pqpddq_poe - p2tau_pqpddq_rst)

norm3d(p2ddq_pqpq_ga - p2ddq_pqpq_rst)
norm3d(p2ddq_pdqpq_ga - p2ddq_pdqpq_rst)
norm3d(p2ddq_pdqpdq_ga - p2ddq_pdqpdq_rst)
norm3d(p2ddq_ptaupq_ga - p2ddq_ptaupq_rst)

norm3d(p2ddq_pqpq_poe - p2ddq_pqpq_rst)
norm3d(p2ddq_pdqpq_poe - p2ddq_pdqpq_rst)
norm3d(p2ddq_pdqpdq_poe - p2ddq_pdqpdq_rst)
norm3d(p2ddq_ptaupq_poe - p2ddq_ptaupq_rst)