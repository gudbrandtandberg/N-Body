%% main.m
% Main plotting program for N-Body Simulations 

function [] = main()

% Load data
positions = load('../output/13_body_trajectories_100_0.0_0.dat');
energies = load('../output/13_body_energy_100_0.0_0.dat');
dt = 0.01;

% Animation of system

%animate(positions, 1:13, 'Solar', 1, 5, 0.1, dt);
%animate(positions, [4 7], 'Earth-Moon', 4, 0.01, 0.001, dt);
%animate(positions, [5 8 9], 'Martian', 5, 0.01, 0.001, dt);
%animate(positions, [6 10 11 12 13], 'Giovian', 6, 0.02, 0.001, dt);

%static_plot(positions);

plot_cluster();

end

% Static plot of entire system
function [] = static_plot(positions)
set(gca, 'FontSize', 16);
set(gcf, 'Color', 'white');
AX = 5;
N = size(positions, 2)/3;
hold on
for i = 1:N
    plot(positions(:, 3*i-2), positions(:, 3*i-1));
end
title('Inner solar system');
xlabel('x [AU]');
ylabel('y [AU]');

axis([-AX AX -AX AX]);
end

% Animates the solar system or a part of the solar system
function [] = animate(positions, I, system, central, ax_size, im_size, dt) 

bodyimg = imread('../images/redbody.png');
n = size(positions, 1); % # of timesteps
set(gca, 'FontSize', 16);
set(gcf, 'Color', 'white');

for i = 1:n
    plot(positions(1:i, I(1)), positions(1:i, I(2)), 'k');
    hold on;
    for j = I
            plot(positions(1:i, 3*j-2), positions(1:i, 3*j-1));
    end
    
    for j = I
        image([positions(i, j*3-2)-im_size positions(i, j*3-2)+im_size],...
            [positions(i, j*3-1)+im_size positions(i, j*3-1)-im_size], bodyimg);
    end
    
    hold off;
    title(sprintf('The %s System %.1f weeks after 1.10.2014', system, dt*i));   
    xlabel('Distance [AU]');
    ylabel('Distance [AU]');
    axis([positions(i, 3*central-2)-ax_size positions(i, 3*central-2)+ax_size...
        positions(i, 3*central-1)-ax_size positions(i, 3*central-1)+ax_size]); 
    pause(0.01);
    hold off;
end
end

function [] = plot_cluster()

positions = load('../output/100_body_trajectories_10_0.1_0.dat');
N = size(positions, 1)/3;
n = size(positions, 2);

for i = 1:n
    hold on;
    
    for j = 1:N
        plot3(positions(1:i, 3*j-2), positions(1:i, 3*j-1), positions(1:i, 3*j));
        
    end
    
    hold off;
    pause(0.001);
    axis([-20 20 -20 20 -20 20]);
end

end


