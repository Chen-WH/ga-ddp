function T = POE_motor(xi, q, type)
% L 是关节轴旋量
% q 是关节运动量
% type 是关节类型:
% 0 - rigidbody motion
% 1 - revolution
% 2 - translation
switch type
    case 0
        Omega = [0, -xi(3),  xi(2);
             xi(3),      0, -xi(1);
            -xi(2),  xi(1),     0];
        R = eye(3) + sin(q)*Omega + (1 - cos(q))*Omega^2;
        A = eye(3) + Omega/q*(1 - cos(q)) + Omega^2/q*(q - sin(q));
        T = [R, A*xi(4:6); 0, 0, 0, 1];
    case 1 % 27M 12A 1sin 1cos
        Omega = [0, -xi(3),  xi(2);
             xi(3),      0, -xi(1);
            -xi(2),  xi(1),     0];
        T = [eye(3) + sin(q)*Omega + (1 - cos(q))*Omega^2, zeros(3, 1); 0, 0, 0, 1];
    case 2 % 3M
        T = [eye(3), q*xi(4:6); 0, 0, 0, 1];
end
end