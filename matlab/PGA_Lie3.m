function dQ = PGA_Lie3(Q, B)
% 对 trivector ["e123", "e023", "e031", "e012"] 和 bivector ["e01", "e02", "e03", "e23", "e31", "e12"] 做李括号运算
% 输出的量是一个 trivector ["e123", "e023", "e031", "e012"]  dQ = -(QB - BQ)/2
% 6M 6A
dQ = [0;
    B(1) - B(5)*Q(4) + B(6)*Q(3);
    B(2) + B(4)*Q(4) - B(6)*Q(2);
    B(3) - B(4)*Q(3) + B(5)*Q(2)];
end