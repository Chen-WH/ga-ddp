function X = CGA_join(type, x)
% 求解多个几何特征交(wedge)的几何元素
% 输入和输出的CGA几何特征均已单位化
% type 10 Sphere ["ei", "eo", "e1", "e2", "e3"]
% type 11 Plane ["ei", "eo", "e1", "e2", "e3"]
% type 20 Circle ["eio", "ei1", "ei2", "ei3", "eo1", "eo2", "eo3", "e23", "e31", "e12"]
% type 21 Line ["eio", "ei1", "ei2", "ei3", "eo1", "eo2", "eo3", "e23", "e31", "e12"]
% type 30 Ponit pair ["eio1", "eio2", "eio3", "ei23", "ei31", "ei12", "eo23", "eo31", "eo12", "e123"]
% type 31 Flat point ["eio1", "eio2", "eio3", "ei23", "ei31", "ei12", "eo23", "eo31", "eo12", "e123"]

switch type
    case 10
        X = [];
    case 11
        X = [];
    case 20
        X = [];
    case 21
        X = [];
    case 30
        X = [];
    case 31
        X = [];
end
end