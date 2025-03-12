function omega = SO3_log(R)
% Log function on SO(3)
% 输入为李群SO(3)，输出为其李代数so(3)
theta = acos((trace(R) - 1)/2);
if abs(theta) < 1e-9
    omega = [0; 0; 0];
else
    Omega = theta*(R - R.')/(2*sin(theta));
    omega = [Omega(3, 2); Omega(1, 3); Omega(2, 1)];
end
end