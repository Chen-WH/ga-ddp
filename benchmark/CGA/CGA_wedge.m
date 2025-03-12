function res = CGA_wedge(l, r, varargin)
% CGA wedge product
%
% basis vector = ["1", "ei", "eo", "e1", "e2", "e3",
% "eio", "ei1", "ei2", "ei3", "eo1", "eo2", "eo3", "e23", "e31", "e12",
% "eio1", "eio2", "eio3", "ei23", "ei31", "ei12", "eo23", "eo31", "eo12", "e123",
% "eio23", "eio31", "eio12", "ei123", "eo123", "eio123"]

res = [l(1)*r(1); % 1
    l(1)*r(2) + l(2)*r(1); % ei
    l(1)*r(3) + l(3)*r(1); % eo
    l(1)*r(4) + l(4)*r(1); % e1
    l(1)*r(5) + l(5)*r(1); % e2
    l(1)*r(6) + l(6)*r(1); % e3

    l(2)*r(3) - l(3)*r(2) + l(1)*r(7) + l(7)*r(1); % eio
    l(2)*r(4) - l(4)*r(2) + l(1)*r(8) + l(8)*r(1); % ei1
    l(2)*r(5) - l(5)*r(2) + l(1)*r(9) + l(9)*r(1); % ei2
    l(2)*r(6) - l(6)*r(2) + l(1)*r(10) + l(10)*r(1); % ei3
    l(3)*r(4) - l(4)*r(3) + l(1)*r(11) + l(11)*r(1); % eo1
    l(3)*r(5) - l(5)*r(3) + l(1)*r(12) + l(12)*r(1); % eo2
    l(3)*r(6) - l(6)*r(3) + l(1)*r(13) + l(13)*r(1); % eo3
    l(5)*r(6) - l(6)*r(5) + l(1)*r(14) + l(14)*r(1); % e23
    l(6)*r(4) - l(4)*r(6) + l(1)*r(15) + l(15)*r(1); % e31
    l(4)*r(5) - l(5)*r(4) + l(1)*r(16) + l(16)*r(1); % e12

    l(4)*r(7) - l(3)*r(8) + l(7)*r(4) - l(8)*r(3) + l(2)*r(11) + l(11)*r(2) + l(1)*r(17) + l(17)*r(1); % eio1
    l(5)*r(7) - l(3)*r(9) + l(7)*r(5) - l(9)*r(3) + l(2)*r(12) + l(12)*r(2) + l(1)*r(18) + l(18)*r(1); % eio2
    l(6)*r(7) - l(3)*r(10) + l(7)*r(6) - l(10)*r(3) + l(2)*r(13) + l(13)*r(2) + l(1)*r(19) + l(19)*r(1); % eio3
    l(6)*r(9) - l(5)*r(10) + l(9)*r(6) - l(10)*r(5) + l(2)*r(14) + l(14)*r(2) + l(1)*r(20) + l(20)*r(1); % ei23
    l(4)*r(10) - l(6)*r(8) - l(8)*r(6) + l(10)*r(4) + l(2)*r(15) + l(15)*r(2) + l(1)*r(21) + l(21)*r(1); % ei31
    l(5)*r(8) - l(4)*r(9) + l(8)*r(5) - l(9)*r(4) + l(2)*r(16) + l(16)*r(2) + l(1)*r(22) + l(22)*r(1); % ei12
    l(3)*r(14) + l(14)*r(3) - l(5)*r(13) + l(6)*r(12) + l(12)*r(6) - l(13)*r(5) + l(1)*r(23) + l(23)*r(1); % eo23
    l(4)*r(13) - l(6)*r(11) - l(11)*r(6) + l(13)*r(4) + l(3)*r(15) + l(15)*r(3) + l(1)*r(24) + l(24)*r(1); % eo31
    l(5)*r(11) - l(4)*r(12) + l(11)*r(5) - l(12)*r(4) + l(3)*r(16) + l(16)*r(3) + l(1)*r(25) + l(25)*r(1); % eo12
    l(4)*r(14) + l(14)*r(4) + l(5)*r(15) + l(15)*r(5) + l(6)*r(16) + l(16)*r(6) + l(1)*r(26) + l(26)*r(1); % e123

    l(7)*r(14) + l(14)*r(7) - l(9)*r(13) + l(10)*r(12) + l(12)*r(10) - l(13)*r(9) - l(3)*r(20) + l(20)*r(3) + l(5)*r(19) - l(6)*r(18) + l(18)*r(6) - l(19)*r(5) + l(2)*r(23) - l(23)*r(2) + l(1)*r(27) + l(27)*r(1); % eio23
    l(8)*r(13) - l(10)*r(11) - l(11)*r(10) + l(13)*r(8) + l(7)*r(15) + l(15)*r(7) - l(4)*r(19) + l(6)*r(17) - l(17)*r(6) + l(19)*r(4) - l(3)*r(21) + l(21)*r(3) + l(2)*r(24) - l(24)*r(2) + l(1)*r(28) + l(28)*r(1); % eio31
    l(9)*r(11) - l(8)*r(12) + l(11)*r(9) - l(12)*r(8) + l(4)*r(18) - l(5)*r(17) + l(17)*r(5) - l(18)*r(4) + l(7)*r(16) + l(16)*r(7) - l(3)*r(22) + l(22)*r(3) + l(2)*r(25) - l(25)*r(2) + l(1)*r(29) + l(29)*r(1); % eio12
    l(8)*r(14) + l(14)*r(8) - l(4)*r(20) + l(9)*r(15) + l(15)*r(9) + l(20)*r(4) - l(5)*r(21) + l(10)*r(16) + l(16)*r(10) + l(21)*r(5) + l(2)*r(26) - l(6)*r(22) + l(22)*r(6) - l(26)*r(2) + l(1)*r(30) + l(30)*r(1); % ei123
    l(11)*r(14) + l(14)*r(11) - l(4)*r(23) + l(12)*r(15) + l(15)*r(12) + l(23)*r(4) + l(3)*r(26) - l(5)*r(24) + l(13)*r(16) + l(16)*r(13) + l(24)*r(5) - l(26)*r(3) - l(6)*r(25) + l(25)*r(6) + l(1)*r(31) + l(31)*r(1); % eo123
    l(4)*r(27) - l(8)*r(23) + l(11)*r(20) + l(14)*r(17) + l(17)*r(14) + l(20)*r(11) - l(23)*r(8) + l(27)*r(4) + l(1)*r(32) + l(2)*r(31) - l(3)*r(30) + l(5)*r(28) + l(7)*r(26) - l(9)*r(24) + l(12)*r(21) + l(15)*r(18) + l(18)*r(15) + l(21)*r(12) - l(24)*r(9) + l(26)*r(7) + l(28)*r(5) - l(30)*r(3) + l(31)*r(2) + l(32)*r(1) + l(6)*r(29) - l(10)*r(25) + l(13)*r(22) + l(16)*r(19) + l(19)*r(16) + l(22)*r(13) - l(25)*r(10) + l(29)*r(6)]; % eio123

if ~isempty(varargin)
    res = CGA_wedge(res, varargin{1}, varargin{2:end});
end

end