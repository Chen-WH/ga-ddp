function res = GA_metric(x)
% 李代数上的度量映射
% res = x \wedge
% basis vector = ["1", "e0", "e1", "e2", "e3", "e01", "e02", "e03", "e23", "e31", "e12", "e123", "e023", "e031", "e012", "e0123"]
res = [x(4), x(5), x(6), x(1), x(2), x(3)];
end