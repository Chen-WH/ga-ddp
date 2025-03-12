function X = CGA_element(type, x)
% basis vector = ["1", "ei", "eo", "e1", "e2", "e3",
% "eio", "ei1", "ei2", "ei3", "eo1", "eo2", "eo3", "e23", "e31", "e12",
% "eio1", "eio2", "eio3", "ei23", "ei31", "ei12", "eo23", "eo31", "eo12", "e123",
% "eio23", "eio31", "eio12", "ei123", "eo123", "eio123"]
%
% type 10 Sphere = [sphere center x; radius rho]
% type 11 Plane = [normal vector n; distance d]
% type 12 Point = [vector x]
% type 20 Circle = [circle center x; radius rho]
% type 21 Line = []
% type 30 Ponit pair = []
% type 31 Flat point = []
% type 40 Point = [vector x]

switch type
    case 10
        X = [(x(1)*x(1) + x(2)*x(2) + x(3)*x(3) - x(4)*x(4))/2; 1; x(1); x(2); x(3)]; % vector = ["ei", "eo", "e1", "e2", "e3"]
    case 11
        X = []; % vector = ["ei", "eo", "e1", "e2", "e3"]
    case 12
        X = [(x(1)*x(1) + x(2)*x(2) + x(3)*x(3))/2; 1; x(1); x(2); x(3)]; % vector = ["ei", "eo", "e1", "e2", "e3"]
    case 20
        X = []; % bivector = ["ei", "eo", "e1", "e2", "e3"]
    case 21
        X = []; % bivector = ["ei", "eo", "e1", "e2", "e3"]
    case 30
        X = []; % trivector = ["ei", "eo", "e1", "e2", "e3"]
    case 31
        X = []; % trivector = ["ei", "eo", "e1", "e2", "e3"]
    case 40
        X = []; % 4-vector = ["eio23", "eio31", "eio12", "ei123", "eo123"]
end
end