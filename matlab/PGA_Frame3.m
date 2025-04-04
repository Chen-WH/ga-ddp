function Qp = PGA_Frame3(Q, M)
% basis vector = ["1", "e0", "e1", "e2", "e3", "e01", "e02", "e03", "e23", "e31", "e12", "e123", "e023", "e031", "e012", "e0123"]
% calculation of M*Q*iM, a rigid motion of unit point Q by M
% 假设Q已单位化
% 58M 36A
Qp = [M(1)^2 + M(5)^2 + M(6)^2 + M(7)^2;
    Q(2)*(M(1)^2 + M(5)^2 - M(6)^2 - M(7)^2) + 2*(-Q(4)*M(1)*M(6) + Q(3)*M(1)*M(7) + M(2)*M(1) + Q(3)*M(5)*M(6) + M(8)*M(5) + Q(4)*M(5)*M(7) - M(4)*M(6) + M(3)*M(7));
    Q(3)*(M(1)^2 - M(5)^2 + M(6)^2 - M(7)^2) + 2*( Q(4)*M(1)*M(5) - Q(2)*M(1)*M(7) + M(3)*M(1) + Q(2)*M(5)*M(6) + M(4)*M(5) + Q(4)*M(6)*M(7) + M(8)*M(6) - M(2)*M(7));
    Q(4)*(M(1)^2 - M(5)^2 - M(6)^2 + M(7)^2) + 2*(-Q(3)*M(1)*M(5) + Q(2)*M(1)*M(6) + M(4)*M(1) + Q(2)*M(5)*M(7) - M(3)*M(5) + Q(3)*M(6)*M(7) + M(2)*M(6) + M(8)*M(7))];
end