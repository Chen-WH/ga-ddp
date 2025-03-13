clc; clear; close all;
rbt_ga = GA_panda();
rbt_rst = importrobot(rbt_ga.urdf, DataFormat='column');
q = rand(rbt_ga.n, 1)*pi - pi/2;
% q = [1.68084; 0.591165; 2.6006; 0.301313; -1.35749; 0.689955];
T = getTransform(rbt_rst, q, rbt_rst.BodyNames{end})

% fkine
M = GA_fkine(rbt_ga, q, rbt_ga.n+1);
O = PGA_Frame3([1; 0; 0; 0], M); O = -O(2:4);
X = PGA_Frame3([1; -1; 0; 0], M); X = -X(2:4) - O;
Y = PGA_Frame3([1; 0; -1; 0], M); Y = -Y(2:4) - O;
Z = PGA_Frame3([1; 0; 0; -1], M); Z = -Z(2:4) - O;
T_ga = [X, Y, Z, O; 0, 0, 0, 1]
norm(T_ga/T)

% ikine
q0 = rand(rbt_ga.n, 1)*pi - pi/2;
q_ga = GA_ikine(rbt_ga, M, q0);
norm(GA_prodM(GA_fkine(rbt_ga, q_ga, rbt_ga.n+1), GA_Rev(M)))