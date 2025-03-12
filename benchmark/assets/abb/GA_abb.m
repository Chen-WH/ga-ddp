function robot = GA_abb()
% The model of the robot
robot.n = 6;
robot.g = [0; 0; -9.81];
robot.urdf = "abb.urdf";

% The motor between adjacent joint frames generated by MDH parameter
robot.Mj(:, 1) = GA_DHmotor([    0;     0; 0.495; 0], 1);
robot.Mj(:, 1) = GA_prodM(GA_exp([-1.3208; -0.9529; -0.7173; 0; 0; 0]), GA_exp([0; 0; 0; 0; 0; pi/2]), robot.Mj(:, 1));
robot.Mj(:, 2) = GA_DHmotor([0.175; -pi/2;     0; -pi/2], 1);
robot.Mj(:, 3) = GA_DHmotor([  0.9;     0;     0; 0], 1);
robot.Mj(:, 4) = GA_DHmotor([0.175; -pi/2;  0.96; 0], 1);
robot.Mj(:, 5) = GA_DHmotor([    0;  pi/2;     0; 0], 1);
robot.Mj(:, 6) = GA_DHmotor([    0; -pi/2;     0; 0], 1);
robot.Mj(:, 7) = GA_DHmotor([    0;     0; 0.135; 0], 1); % end effector is frame {i+1}
% robot.Mj(:, 7) = GA_DHmotor([    0;     0;  0.16; 0], 1); % with target balls

% The initial position of the joint axis
for i = 1:robot.n
    robot.Lj(:, i) = [0; 0; 0; 0; 0; 1]; % in the joint frame
    robot.type(i) = 1;
end