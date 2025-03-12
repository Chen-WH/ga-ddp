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
[ptau_pq_ga, ptau_pdq_ga, ptau_pddq_ga] = GA_idyn_fo(rbt_ga, q, dq, ddq);
[pddq_pq_ga, pddq_pdq_ga, pddq_ptau_ga] = GA_fdyn_fo(rbt_ga, q, dq, tau);

%% POE
[ptau_pq_poe, ptau_pdq_poe, ptau_pddq_poe] = POE_idyn_fo(rbt_poe, q, dq, ddq);
[pddq_pq_poe, pddq_pdq_poe, pddq_ptau_poe] = POE_fdyn_fo(rbt_poe, q, dq, tau);

%% Finite difference
ptau_pq_rst = zeros(rbt_ga.n, rbt_ga.n);
ptau_pdq_rst = zeros(rbt_ga.n, rbt_ga.n);
ptau_pddq_rst = zeros(rbt_ga.n, rbt_ga.n);
pddq_pq_rst = zeros(rbt_ga.n, rbt_ga.n);
pddq_pdq_rst = zeros(rbt_ga.n, rbt_ga.n);
pddq_ptau_rst = zeros(rbt_ga.n, rbt_ga.n);

step = 1e-7;
for i = 1:rbt_ga.n
    q_tmp = q;
    dq_tmp = dq;
    ddq_tmp = ddq;
    tau_tmp = tau;
    q_tmp(i) = q_tmp(i) + step;
    dq_tmp(i) = dq_tmp(i) + step;
    ddq_tmp(i) = ddq_tmp(i) + step;
    tau_tmp(i) = tau_tmp(i) + step;
    % ID-FO
    ptau_pq_rst(:, i) = (inverseDynamics(rbt_rst, q_tmp, dq, ddq) - tau)/step;
    ptau_pdq_rst(:, i) = (inverseDynamics(rbt_rst, q, dq_tmp, ddq) - tau)/step;
    ptau_pddq_rst(:, i) = (inverseDynamics(rbt_rst, q, dq, ddq_tmp) - tau)/step;
    % FD-FO
    pddq_pq_rst(:, i) = (forwardDynamics(rbt_rst, q_tmp, dq, tau) - ddq)/step;
    pddq_pdq_rst(:, i) = (forwardDynamics(rbt_rst, q, dq_tmp, tau) - ddq)/step;
    pddq_ptau_rst(:, i) = (forwardDynamics(rbt_rst, q, dq, tau_tmp) - ddq)/step;
end

%% Result Validation
norm(ptau_pq_ga - ptau_pq_rst)
norm(ptau_pdq_ga - ptau_pdq_rst)
norm(ptau_pddq_ga - ptau_pddq_rst)

norm(ptau_pq_poe - ptau_pq_rst)
norm(ptau_pdq_poe - ptau_pdq_rst)
norm(ptau_pddq_poe - ptau_pddq_rst)

norm(pddq_pq_ga - pddq_pq_rst)
norm(pddq_pdq_ga - pddq_pdq_rst)
norm(pddq_ptau_ga - pddq_ptau_rst)

norm(pddq_pq_poe - pddq_pq_rst)
norm(pddq_pdq_poe - pddq_pdq_rst)
norm(pddq_ptau_poe - pddq_ptau_rst)