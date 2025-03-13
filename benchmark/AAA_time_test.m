clc; clear; close all;
rbt_ga = GA_panda();
rbt_sva = SVA_panda();
rbt_poe = POE_panda();
rbt_rst = importrobot(rbt_ga.urdf, DataFormat='column');
rbt_rst.Gravity = rbt_ga.g;
n = 10000;
q = rand(rbt_ga.n, n)*2*pi - pi;
dq = rand(rbt_ga.n, n)*10 - 5;
ddq = rand(rbt_ga.n, n)*10 - 5;
tau = zeros(rbt_ga.n, n);
for i = 1:n
    tau(:, i) = inverseDynamics(rbt_rst, q(:, i), dq(:, i), ddq(:, i));
end

%% GA
tic
for i = 1:n
    GA_fdyn(rbt_ga, q(:, i), dq(:, i), tau(:, i));
end
Time(1, 1) = toc/n;

tic
for i = 1:n
    GA_fdyn_fo(rbt_ga, q(:, i), dq(:, i), tau(:, i));
end
Time(2, 1) = toc/n;

tic
for i = 1:n
    GA_fdyn_so(rbt_ga, q(:, i), dq(:, i), tau(:, i));
end
Time(3, 1) = toc/n;

%% SVA
tic
for i = 1:n
    FDab(rbt_sva, q(:, i), dq(:, i), tau(:, i));
end
Time(1, 2) = toc/n;

tic
for i = 1:n
    modFD_derivatives(rbt_sva, q(:, i), dq(:, i), tau(:, i), rand(rbt_ga.n, 1));
end
Time(2, 2) = toc/n;

tic
for i = 1:n
    modFD_second_derivatives(rbt_sva, q(:, i), dq(:, i), tau(:, i), rand(rbt_ga.n, 1));
end
Time(3, 2) = toc/n;

%% POE
tic
for i = 1:n
    POE_fdyn(rbt_poe, q(:, i), dq(:, i), tau(:, i));
end
Time(1, 3) = toc/n;

tic
for i = 1:n
    POE_fdyn_fo(rbt_poe, q(:, i), dq(:, i), tau(:, i));
end
Time(2, 3) = toc/n;

tic
for i = 1:n
    POE_fdyn_so(rbt_poe, q(:, i), dq(:, i), tau(:, i));
end
Time(3, 3) = toc/n;

%% Figure
Time = 1000*Time;

bar(["Forward Dynamics", "First-order Derivs", "Second-order Derivs"], Time);
ylabel("Time(ms)");
legend("GA", "SVA", "POE", "Location","best");
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',11);
saveas(gcf, "../images/benchmark.png");