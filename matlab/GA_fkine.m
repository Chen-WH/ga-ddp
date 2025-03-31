function M = GA_fkine(rbt, q, sourcebody, targetbody)
if nargin == 3
    targetbody = 0;
end
n = rbt.n;
M0 = rbt.Mj;
L = rbt.Lj;
type = rbt.type;
if sourcebody == n+1
    M = M0(:, n+1);
else
    M = [1; zeros(7, 1)];
end
for i = min(sourcebody, n):-1:targetbody+1
    M = GA_prodM(M0(:, i), GA_motor(L(:, i), q(i), type(i)), M);
end
end