% Static plot of entire system
function [] = static_plot(positions, tit)
set(gca, 'FontSize', 16);
set(gcf, 'Color', 'white');
AX = 6;
N = size(positions, 2)/3;
hold on
axis equal
for i = 1:N
    plot(positions(:, 3*i-2), positions(:, 3*i-1));
end
title(tit);
xlabel('x [AU]');
ylabel('y [AU]');

axis([-AX AX -AX AX]);
end
