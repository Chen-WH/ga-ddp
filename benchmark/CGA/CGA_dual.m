function X_dual = CGA_dual(X)
% basis vector = ["1", "ei", "eo", "e1", "e2", "e3",
% "eio", "ei1", "ei2", "ei3", "eo1", "eo2", "eo3", "e23", "e31", "e12",
% "eio1", "eio2", "eio3", "ei23", "ei31", "ei12", "eo23", "eo31", "eo12", "e123",
% "eio23", "eio31", "eio12", "ei123", "eo123", "eio123"]

n = length(X);
switch n
    case 32
        X_dual = [X(32); X(30); -X(31); X(27); X(28); X(29);
            X(26); X(20); X(21); X(22); -X(23); -X(24); -X(25); -X(17); -X(18); -X(19);
            X(14); X(15); X(16); -X(8); -X(9); -X(10); X(11); X(12); X(13); -X(7);
            -X(4); -X(5); -X(6); -X(2); X(3); -X(1)];
end
end