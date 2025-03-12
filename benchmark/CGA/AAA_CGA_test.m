clc;clear;close all;
% basis vector = ["1", "ei", "eo", "e1", "e2", "e3",
% "eio", "ei1", "ei2", "ei3", "eo1", "eo2", "eo3", "e23", "e31", "e12",
% "eio1", "eio2", "eio3", "ei23", "ei31", "ei12", "eo23", "eo31", "eo12", "e123",
% "eio23", "eio31", "eio12", "ei123", "eo123", "eio123"]
% conformal_signature(3, 0);

%% 推导 CGA element
syms x11 x31 x41 x51 x12 x32 x42 x52 x13 x33 x43 x53 x14 x34 x44 x54
einfty = [0; 1; zeros(30, 1)];
X1 = [0; x11; 1; x31; x41; x51; zeros(26, 1)];
X2 = [0; x12; 1; x32; x42; x52; zeros(26, 1)];
X3 = [0; x13; 1; x33; x43; x53; zeros(26, 1)];
X4 = [0; x14; 1; x34; x44; x54; zeros(26, 1)];
combine(CGA_dual(CGA_wedge(X1, X2, X3, einfty)))