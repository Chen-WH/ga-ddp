function M = GA_motor(L, q, type)
% M = exp(-qL/2)
% L 是关节轴 bivector, |L(4:6)|=1 已单位化
% q 是关节量 scalar
% type 是关节类型:
% 0 - screw, q = [d; theta]
% 1 - revolution, q = theta
% 2 - translation, q = d
q = -0.5*q;
switch type
    case 0 % 13M 3A 1sin 1cos
        c = cos(q(2));
        s = sin(q(2));
        M = [c; s*L(1:3) + c*q(1)*L(4:6); s*L(4:6); s*q(1)];
    case 1 % 7M 1sin 1cos
        M = [cos(q); sin(q)*L; 0];
    case 2 % 4M
        M = [1; q*L(1:3); 0; 0; 0; 0];
end
end