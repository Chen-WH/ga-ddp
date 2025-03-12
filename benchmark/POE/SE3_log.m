function xi = SE3_log(T)
% Log function on SE(3)
% 输入为李群SE(3)，输出为其李代数se(3)
R = T(1:3, 1:3);
p = T(1:3, 4);
theta = acos((trace(R) - 1)/2);
if abs(theta) < 1e-9
    xi = [0; 0; 0; p];
else
    Omega = (R - R.')/(2*sin(theta));
    omega = [Omega(3, 2); Omega(1, 3); Omega(2, 1)];
    invA = eye(3) - theta*Omega/2 + (2*sin(theta) - theta*(1 + cos(theta)))/(2*sin(theta))*Omega^2;
    xi = [theta*omega; invA*p];
end
end