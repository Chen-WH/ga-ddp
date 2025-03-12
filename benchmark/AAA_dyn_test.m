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
tau_ga = GA_idyn(rbt_ga, q, dq, ddq);
ddq_ga = GA_fdyn(rbt_ga, q, dq, tau);

%% POE
tau_poe = POE_idyn(rbt_poe, q, dq, ddq);
ddq_poe = POE_fdyn(rbt_poe, q, dq, tau);

%% Result Validation
norm(tau - tau_poe)
norm(tau - tau_ga)
norm(ddq - ddq_poe)
norm(ddq - ddq_ga)