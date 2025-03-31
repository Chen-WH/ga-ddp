function dM = GA_prodMB(M, Bb)
% 对 motor 和 bivector ["e01", "e02", "e03", "e23", "e31", "e12"] 做几何积运算
% 输出的量是motor的变化率
% dM = M*Bb (36M 28A)
dM = [- Bb(4)*M(5) - Bb(5)*M(6) - Bb(6)*M(7);
    Bb(1)*M(1) + Bb(2)*M(7) - Bb(3)*M(6) + Bb(5)*M(4) - Bb(6)*M(3) - Bb(4)*M(8);
    Bb(2)*M(1) - Bb(1)*M(7) + Bb(3)*M(5) - Bb(4)*M(4) + Bb(6)*M(2) - Bb(5)*M(8);
    Bb(3)*M(1) + Bb(1)*M(6) - Bb(2)*M(5) + Bb(4)*M(3) - Bb(5)*M(2) - Bb(6)*M(8);
    Bb(4)*M(1) + Bb(5)*M(7) - Bb(6)*M(6);
    Bb(5)*M(1) - Bb(4)*M(7) + Bb(6)*M(5);
    Bb(6)*M(1) + Bb(4)*M(6) - Bb(5)*M(5);
    Bb(1)*M(5) + Bb(4)*M(2) + Bb(2)*M(6) + Bb(5)*M(3) + Bb(3)*M(7) + Bb(6)*M(4)];
end