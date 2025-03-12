function M = GA_DHmotor(param, type)
% DH 参数
% type 是 DH 存储类型
% 0 -  DH, param = [d, theta, a, alpha], M = exp(d*e03 + theta*e12)*exp(a*e01 + alpha*e23)
% 1 - MDH, param = [a, alpha, d, theta], M = exp(a*e01 + alpha*e23)*exp(d*e03 + theta*e12)
param = 0.5*param;
switch type
    case 0
        d = param(1); theta = param(2); a = param(3); alpha = param(4);
        M = [cos(alpha)*cos(theta);
            d*sin(alpha)*sin(theta) - a*cos(alpha)*cos(theta);
           -a*cos(alpha)*sin(theta) - d*sin(alpha)*cos(theta);
            a*sin(alpha)*sin(theta) - d*cos(alpha)*cos(theta);
           -sin(alpha)*cos(theta);
           -sin(alpha)*sin(theta);
           -cos(alpha)*sin(theta);
            a*sin(alpha)*cos(theta) + d*cos(alpha)*sin(theta)];
    case 1
        a = param(1); alpha = param(2); d = param(3); theta = param(4);
        M = [cos(alpha)*cos(theta);
            d*sin(alpha)*sin(theta) - a*cos(alpha)*cos(theta);
            a*cos(alpha)*sin(theta) + d*sin(alpha)*cos(theta);
            a*sin(alpha)*sin(theta) - d*cos(alpha)*cos(theta);
            -sin(alpha)*cos(theta);
            sin(alpha)*sin(theta);
            -cos(alpha)*sin(theta);
            a*sin(alpha)*cos(theta) + d*cos(alpha)*sin(theta)];
end
end