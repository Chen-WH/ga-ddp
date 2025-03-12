function X = PGA_element(type, coord)
% basis vector = ["1", "e0", "e1", "e2", "e3", "e01", "e02", "e03", "e23", "e31", "e12", "e123", "e023", "e031", "e012", "e0123"]
%
% type 10 Plane = [normal vector n; distance delta]
% type 20 Line = [vector q; vector u]
% type 21 Hyperline = [normal vector n1; distance delta1 ; normal vector n2; distance delta2]
% type 30 Point = [vector q]
% type 31 Mass Point = [vector q; mass m]
switch type
    case 10
        X = [-coord(4); coord(1); coord(2); coord(3)]; % vector = ["e0", "e1", "e2", "e3"]
    case 20
        X = [coord(2)*coord(6) - coord(3)*coord(5); coord(3)*coord(4) - coord(1)*coord(6); coord(1)*coord(5) - coord(2)*coord(4);
            coord(4); coord(5); coord(6)]; % bivector = ["e01", "e02", "e03", "e23", "e31", "e12"]
    case 21
        X = [coord(1)*coord(8) - coord(4)*coord(5); coord(2)*coord(8) - coord(4)*coord(6); coord(3)*coord(8) - coord(4)*coord(7);
            coord(2)*coord(7) - coord(3)*coord(6); coord(3)*coord(5) - coord(1)*coord(7); coord(1)*coord(6) - coord(2)*coord(5)]; % bivector = ["e01", "e02", "e03", "e23", "e31", "e12"]
    case 30
        X = [1; -coord(1); -coord(2); -coord(3)]; % trivector = ["e123", "e023", "e031", "e012"]
    case 31
        X = [coord(4); -coord(4)*coord(1); -coord(4)*coord(2); -coord(4)*coord(3)]; % trivector = ["e123", "e023", "e031", "e012"]
end
end