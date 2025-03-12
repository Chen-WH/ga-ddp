function R = SO3_exp(omega)
% Exp function on SO(3)
% 输入为李代数so(3)，输出为其指数映射李群SO(3)
theta = norm(omega);
omega = omega/theta;
Omega = [0, -omega(3),  omega(2);
  omega(3),         0, -omega(1);
 -omega(2),  omega(1),        0];
R = eye(3) + sin(theta)*Omega + (1 - cos(theta))*Omega^2;
end