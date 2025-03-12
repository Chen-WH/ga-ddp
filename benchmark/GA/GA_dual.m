function res = GA_dual(x)
% basis vector = ["1", "e0", "e1", "e2", "e3", "e01", "e02", "e03", "e23", "e31", "e12", "e123", "e023", "e031", "e012", "e0123"]
n = length(x);
switch n
    case 6
        res = -[x(4:6); x(1:3)];
    case 8
        res = [x(8); -x(5:7); -x(2:4); x(1)];
end
end