function F = GA_forque(q, f)
% 输入为单位质量受力点 vector q 和 外力 vector f ["e1", "e2", "e3"]
% 输出为forque bivector [torque; force] ["e01", "e02", "e03", "e23", "e31", "e12"]
% F = (O + qI)v(fI) (6M 3A)
F = [q(2)*f(3) - q(3)*f(2);
    q(3)*f(1) - q(1)*f(3);
    q(1)*f(2) - q(2)*f(1);
    f(1);
    f(2);
    f(3)];
end