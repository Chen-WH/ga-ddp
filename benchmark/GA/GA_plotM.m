function GA_plotM(M)

origin = PGA_Frame3([1; 0; 0; 0], M);
x_axis = PGA_Frame3([1; 0.1; 0; 0], M) - origin;
y_axis = PGA_Frame3([1; 0; 0.1; 0], M) - origin;
z_axis = PGA_Frame3([1; 0; 0; 0.1], M) - origin;

quiver3(-origin(2), -origin(3), -origin(4), x_axis(2), x_axis(3), x_axis(4), 'r', 'LineWidth', 1.5); % X轴
quiver3(-origin(2), -origin(3), -origin(4), y_axis(2), y_axis(3), y_axis(4), 'g', 'LineWidth', 1.5); % Y轴
quiver3(-origin(2), -origin(3), -origin(4), z_axis(2), z_axis(3), z_axis(4), 'b', 'LineWidth', 1.5); % Z轴

end