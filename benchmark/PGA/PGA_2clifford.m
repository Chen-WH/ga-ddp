function Y = PGA_2clifford(X)
% 将数组形式的PGA元素转换为matlab clifford形式
Y = X(1) + X(2)*e4 + X(3)*e1 + X(4)*e2 + X(5)*e3 - X(6)*e14 -X(7)*e24 - X(8)*e34 + X(9)*e23 - X(10)*e13 + X(11)*e12 + X(12)*e123 + X(13)*e234 - X(14)*e134 + X(15)*e124 - X(16)*e1234;
end