function T = POE_DHmotor(param, type)
% DH 参数
% type 是 DH 存储类型
% 0 -  DH, param = [d, theta, a, alpha], M = exp(alpha*e23 + a*e01)*exp(theta*e12 + d*e03)
% 1 - MDH, param = [a, alpha, d, theta], M = exp(theta*e12 + d*e03)*exp(alpha*e23 + a*e01)
switch type
    case 0
        d = param(1); theta = param(2); a = param(3); alpha = param(4);
        T = [cos(theta), -sin(theta)*cos(alpha),  sin(theta)*sin(alpha), a*cos(theta);
             sin(theta),  cos(theta)*cos(alpha), -cos(theta)*sin(alpha), a*sin(theta);
                      0,             sin(alpha),             cos(alpha),            d;
                      0,                      0,                      0,           1];
    case 1
        a = param(1); alpha = param(2); d = param(3); theta = param(4);
        T = [cos(theta),           -sin(theta),           0,             a;
  sin(theta)*cos(alpha), cos(theta)*cos(alpha), -sin(alpha), -sin(alpha)*d;
  sin(theta)*sin(alpha), cos(theta)*sin(alpha),  cos(alpha),  cos(alpha)*d;
                      0,                     0,           0,            1];
end
end