%% main.m
% Main plotting program for N-Body Simulations 

function [] = main()

% Load data

positions = load('../output/250_body_trajectories_5_0.010_1_0_0_0.10.dat');
energies = load('../output/250_body_energy_5_0.010_1_0_0_0.10.dat');
virialenergy = load('../output/250_body_virialenergy_5_0.010_1_0_0_0.10.dat');
numbound = load('../output/250_body_bound_5_0.010_1_0_0_0.10.dat');

plot(virialenergy);
figure();
plot(numbound);

%static_plot(positions, 'Adaptive Verlet, dt_{max} = 1.0');
%figure();
%plot(energies);
%histogram(positions(end, :));

% Vpositions1 = load('../output/6_body_trajectories_619_2.0_0_0.dat');
% Vpositions2 = load('../output/6_body_trajectories_619_1.0_0_0.dat');
% Vpositions3 = load('../output/6_body_trajectories_619_0.1_0_0.dat');
% 
% RK4positions1 = load('../output/6_body_trajectories_619_2.0_0_1.dat');
% RK4positions2 = load('../output/6_body_trajectories_619_1.0_0_1.dat');
% RK4positions3 = load('../output/6_body_trajectories_619_0.1_0_1.dat');
% 
% RK4energies1 = load('../output/6_body_energy_619_2.0_0_1.dat');
% RK4energies2 = load('../output/6_body_energy_619_1.0_0_1.dat');
% RK4energies3 = load('../output/6_body_energy_619_0.1_0_1.dat');
% 
% Venergies1 = load('../output/6_body_energy_619_2.0_0_0.dat');
% Venergies2 = load('../output/6_body_energy_619_1.0_0_0.dat');
% Venergies3 = load('../output/6_body_energy_619_0.1_0_0.dat');
% 
% n1 = size(RK4energies1, 1);
% n2 = size(RK4energies2, 1);
% n3 = size(RK4energies3, 1);

% T = 619;
% t1 = linspace(0, 619, n1);
% t2 = linspace(0, 619, n2);
% t3 = linspace(0, 619, n3);

% figure(1);
% title('Inner solar system trajectories - Verlet method');
% subplot(1, 3, 1);
% static_plot(Vpositions1, 'dt = 2.0');
% subplot(1, 3, 2);
% static_plot(Vpositions2, 'dt = 1.0');
% subplot(1, 3, 3);
% static_plot(Vpositions3, 'dt = 0.1');
% 
% figure(2);
% title('Inner solar system trajectories - RK4 method');
% subplot(1, 3, 1);
% static_plot(RK4positions1, 'dt = 2.0');
% subplot(1, 3, 2);
% static_plot(RK4positions2, 'dt = 1.0');
% subplot(1, 3, 3);
% static_plot(RK4positions3, 'dt = 0.1');
% 
% figure(3);
% 
% plot(t1, Venergies1, t2, Venergies2, 'k', t3, Venergies3, 'r');
% xlabel('Time [weeks]');
% ylabel('Energy');
% legend('dt = 2.0', 'dt = 1.0', 'dt = 0.1');
% title('Energy evolution of inner solar system - Verlet method');
% 
% figure(4);
% 
% plot(t1, RK4energies1, t2, RK4energies2, 'k', t3, RK4energies3, 'r');
% xlabel('Time [weeks]');
% ylabel('Energy');
% legend('dt = 2.0', 'dt = 1.0', 'dt = 0.1');
% title('Energy evolution of inner solar system - RK4 method');
% 
% dt = 0.1;

% Animation of system

%animate(positions, 1:13, 'Solar', 1, 5, 0.05, dt);
%animate(positions, [4 7], 'Earth-Moon', 4, 0.01, 0.001, dt);
%animate(positions, [5 8 9], 'Martian', 5, 0.01, 0.001, dt);
%animate(positions, [6 10 11 12 13], 'Giovian', 6, 0.02, 0.001, dt);


end

% Static plot of entire system
function [] = static_plot(positions, tit)
set(gca, 'FontSize', 16);
set(gcf, 'Color', 'white');
AX = 6;
N = size(positions, 2)/3;
hold on
for i = 1:N
    plot(positions(:, 3*i-2), positions(:, 3*i-1));
end
title(tit);
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

function histogram(positions)
    N = size(positions, 2)/3;
    r = zeros(1, N);
    for i = 1:N
        r(i) = norm([positions(3*i-2) positions(3*i-1) positions(3*i)]);
    end
    
    [dist, centers] = hist(r, 20);
    half_width = (centers(2)-centers(1))/2;
    
    density = zeros(1, 20);
    
    
    for i=1:20
        shell_volume = 4*pi*((centers(i)+half_width)^3 - (centers(i)-half_width)^3)/3;
        density(i) = dist(i)/shell_volume;
    end
    
  
    
    figure();
    n_0 = 200;
    r_0 = 20;
    
    %n_r = n_0./(1+(n./r_0).^4);
 
    plot(density);
    hold on;
    %plot(n_r, 'k');

end


