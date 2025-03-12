function q = GA_ikine(rbt, M, q0)
if nargin == 3
    q = q0;
else
    q = rand(rbt.n, 1);
end
M0 = GA_fkine(rbt, q, rbt.n+1);
e = 1;
iter = 0;
while norm(e) > 1e-5
    e = M0 - M;
    J = GA_jacob_A(rbt, q);
    dq = (J.'*J + 1e-9*eye(rbt.n))\(J.'*e);
    flag = true;
    alpha = 1;
    while flag && alpha > 1e-5
        q_tmp = wrapToPi(q + alpha*dq);
        M_tmp = GA_fkine(rbt, q_tmp, rbt.n+1);
        if norm(M_tmp - M) < norm(e)
            flag = false;
        else
            alpha = 0.618*alpha;
        end
    end
    q = wrapToPi(q + dq);
    M0 = GA_fkine(rbt, q, rbt.n+1);
    if iter > 50
        error("ikine not converge!");
    end
    iter = iter + 1;
end
end