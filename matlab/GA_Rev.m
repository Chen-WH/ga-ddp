function iM = GA_Rev(M)
% basis vector = ["1", "e01", "e02", "e03", "e23", "e31", "e12", "e0123"]
iM = M;
iM(2:7) = -M(2:7);
end