function T = SE3_exp(xi)
% Exp function on SE(3)
% 输入为李代数se(3)，输出为其指数映射李群SE(3)
omega = xi(1:3);
v = xi(4:6);
theta = norm(omega);
if theta < 1e-9
    T = [eye(3), v; 0, 0, 0, 1];
else
    omega = omega/theta;
    Omega = [0, -omega(3),  omega(2);
      omega(3),         0, -omega(1);
     -omega(2),  omega(1),        0];
    R = eye(3) + sin(theta)*Omega + (1 - cos(theta))*Omega^2;
    A = eye(3) + Omega/theta*(1 - cos(theta)) + Omega^2/theta*(theta - sin(theta));
    T = [R, A*v; 0, 0, 0, 1];
end
end