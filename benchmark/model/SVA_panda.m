function robot = SVA_panda()
params = zeros(10, 7);

% mass and inertial parameters of each link in the mass center frame
Nc(:, 1) = [    4.970684;    0.70337; -0.00013900;   0.0067720;    0.70661;    0.019169;  0.0091170];
Nc(:, 2) = [    0.646926;  0.0079620;  -3.9250e-3;  1.0254e-02; 2.8110e-02;  7.0400e-04; 2.5995e-02];
Nc(:, 3) = [    3.228604; 3.7242e-02; -4.7610e-03; -1.1396e-02; 3.6155e-02; -1.2805e-02; 1.0830e-02];
Nc(:, 4) = [    3.587895; 2.5853e-02;  7.7960e-03; -1.3320e-03; 1.9552e-02;  8.6410e-03; 2.8323e-02];
Nc(:, 5) = [    1.225946; 3.5549e-02; -2.1170e-03; -4.0370e-03; 2.9474e-02;  2.2900e-04; 8.6270e-03];
Nc(:, 6) = [    1.666555; 1.9640e-03;  1.0900e-04; -1.1580e-03; 4.3540e-03;  3.4100e-04; 5.4330e-03];
Nc(:, 7) = [ 7.35522e-01; 1.2516e-02; -4.2800e-04; -1.1960e-03; 1.0027e-02; -7.4100e-04; 4.8150e-03];

% Mass center of each link in the joint frame
rc(:, 1) = [    0.003875;    0.002081;    -0.04762];
rc(:, 2) = [   -0.003141;    -0.02872;    0.003495];
rc(:, 3) = [  2.7518e-02;  3.9252e-02; -6.6502e-02];
rc(:, 4) = [  -5.317e-02; 1.04419e-01;  2.7454e-02];
rc(:, 5) = [ -1.1953e-02;  4.1065e-02; -3.8437e-02];
rc(:, 6) = [  6.0149e-02; -1.4117e-02; -1.0517e-02];
rc(:, 7) = [  1.0517e-02;  -4.252e-03;  6.1597e-02];

% calculate  dynamic parameters in the joint frame
for i = 1:7
    I = [Nc(2, i), Nc(3, i), Nc(4, i); Nc(3, i), Nc(5, i), Nc(6, i); Nc(4, i), Nc(6, i), Nc(7, i)];
    Rc = skew(rc(:, i));
    J = I - Nc(1, i)*Rc*Rc;
    params(:, i) = [Nc(1, i); Nc(1, i)*rc(:, i); J(1, 1); J(2, 2); J(3, 3); J(1, 2); J(1, 3); J(2, 3)];
end

robot = Panda_model(params);
robot = postProcessModel(robot);
end