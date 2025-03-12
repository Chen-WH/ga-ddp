function invT = SE3_inv(T)
% Inverse Matrix of Homogeneous Transformation Matrix
% 输入为齐次变换矩阵，输出为其逆矩阵
R_AB = T(1:3,1:3); % 提取旋转矩阵
o_AB = T(1:3,4); % 提取原点坐标
invT = [R_AB.', -R_AB.'*o_AB; 0 0 0 1];
end