% main.m
% Main plotting program for N-Body Simulations 

%% Part 1 - inner solar system tests and trajectories

positions = load('/Users/gudbrand/Documents/C++/build-NBODY-GCC49-Release/output/6_body_trajectories_201_1_0_4_0.00.dat');

static_plot(positions, 'title');

%% Load data
Vpositions1 = load('../output/6_body_trajectories_619_2.00_0_0.dat');
Vpositions2 = load('../output/6_body_trajectories_619_1.00_0_0.dat');
Vpositions3 = load('../output/6_body_trajectories_619_0.05_0_0.dat');

RK4positions1 = load('../output/6_body_trajectories_619_2.00_0_1.dat');
RK4positions2 = load('../output/6_body_trajectories_619_1.00_0_1.dat');
RK4positions3 = load('../output/6_body_trajectories_619_0.05_0_1.dat');

Venergies1 = load('../output/6_body_energy_619_2.00_0_0.dat');
Venergies2 = load('../output/6_body_energy_619_1.00_0_0.dat');
Venergies3 = load('../output/6_body_energy_619_0.05_0_0.dat');

RK4energies1 = load('../output/6_body_energy_619_2.00_0_1.dat');
RK4energies2 = load('../output/6_body_energy_619_1.00_0_1.dat');
RK4energies3 = load('../output/6_body_energy_619_0.05_0_1.dat');

%% Positions
T = 619;

n1 = size(RK4energies1, 1);
n2 = size(RK4energies2, 1);
n3 = size(RK4energies3, 1);

t1 = linspace(0, T, n1);
t2 = linspace(0, T, n2);
t3 = linspace(0, T, n3);

figure();
subplot(1, 3, 1);
static_plot(Vpositions1, 'dt = 2.0');
subplot(1, 3, 2);
static_plot(Vpositions2, {'Inner solar system trajectories - Verlet method', 'dt = 1.0'});
subplot(1, 3, 3);
static_plot(Vpositions3, 'dt = 0.05');

figure();
subplot(1, 3, 1);
static_plot(RK4positions1, 'dt = 2.0');
subplot(1, 3, 2);
static_plot(RK4positions2, {'Inner solar system trajectories - RK4 method', 'dt = 1.0'});
subplot(1, 3, 3);
static_plot(RK4positions3, 'dt = 0.05');

%% Energies
figure();
set(gca, 'FontSize', 16);
set(gcf, 'Color', 'white');
plot(t1, Venergies1, t2, Venergies2, 'k', t3, Venergies3, 'r');
xlabel('Time');
ylabel('Energy');
legend('dt = 2.0', 'dt = 1.0', 'dt = 0.05');
title({'Energy evolution of inner solar system', 'Verlet method'});

figure();
set(gca, 'FontSize', 16);
set(gcf, 'Color', 'white');
plot(t1, RK4energies1, t2, RK4energies2, 'k', t3, RK4energies3, 'r');
xlabel('Time');
ylabel('Energy');
legend('dt = 2.0', 'dt = 1.0', 'dt = 0.05');
title({'Energy evolution of inner solar system', 'RK4 method'});

%% Long simulations

Vpositions2 = load('../output/6_body_trajectories_2000_0.50_0_0.dat');
RK4positions2 = load('../output/6_body_trajectories_2000_0.50_0_1.dat');

figure();
static_plot(Vpositions2, 'V');
figure();
static_plot(RK4positions2, 'RK');

Venergy2 = load('../output/6_body_energy_2000_0.50_0_0.dat');
RK4energy2 = load('../output/6_body_energy_2000_0.50_0_1.dat');
t = linspace(0, 2000, size(Venergy2, 1));

figure();
set(gca, 'FontSize', 16);
set(gcf, 'Color', 'white');
plot(t, Venergy2, t, RK4energy2);
xlabel('Time');
ylabel('Energy');
title({'Long term energy evolution of inner solar system', 'dt = 0.1'});
legend('Verlet', 'RK4');

%% Adaptive solar system (dt = 0.1931, 0.0965, 0.0483)

Apositions1 = load('../output/6_body_trajectories_619_1_0_0_0.00.dat');
Aenergy1 = load('../output/6_body_energy_619_1_0_0_0.00.dat');
plot(Aenergy1);
figure();
static_plot(Apositions1, 'tit');


%% Animations
dt = 0.01;
positions_11 = load('../output/11_body_trajectories_619_0.010_0_0.dat');

figure();
subplot(1, 3, 1);
animate(positions_11, 1:11, 'Inner Solar', 1, 6, 0.05, dt);

subplot(1, 3, 2);
animate(positions_11, [4 7], 'Earth-Moon', 4, 0.01, 0.001, dt);

subplot(1, 3, 3);
animate(positions_11, [6 8 9 10 11], 'Giovian', 6, 0.02, 0.001, dt);


%% Part 3 - star cluster studies

energy1 = load('../output/100_body_energy_6_1_0_1_0.10.dat');
virialenergy1 = load('../output/100_body_virialenergy_6_1_0_1_0.10.dat');
numbound1 = load('../output/100_body_bound_6_1_0_1_0.10.dat')/100;

energy2 = load('../output/250_body_energy_6_1_0_1_0.10.dat');
virialenergy2 = load('../output/250_body_virialenergy_6_1_0_1_0.10.dat');
numbound2 = load('../output/250_body_bound_6_1_0_1_0.10.dat')/250;

energy3 = load('../output/1000_body_energy_6_1_0_1_0.10.dat');
virialenergy3 = load('../output/1000_body_virialenergy_6_1_0_1_0.10.dat');
numbound3 = load('../output/1000_body_bound_6_1_0_1_0.10.dat')/1000;

energy4 = load('../output/100_body_energy_6_1_0_1_0.01.dat');
virialenergy4 = load('../output/100_body_virialenergy_6_1_0_1_0.01.dat');
numbound4 = load('../output/100_body_bound_6_1_0_1_0.01.dat')/100;

energy5 = load('../output/250_body_energy_6_1_0_1_0.01.dat');
virialenergy5 = load('../output/250_body_virialenergy_6_1_0_1_0.01.dat');
numbound5 = load('../output/250_body_bound_6_1_0_1_0.01.dat')/250;

energy6 = load('../output/1000_body_energy_6_1_0_1_0.01.dat');
virialenergy6 = load('../output/1000_body_virialenergy_6_1_0_1_0.01.dat');
numbound6 = load('../output/1000_body_bound_6_1_0_1_0.01.dat')/1000;

T = 6; 

dt1 = 0.0026; n1 = size(numbound1, 1); t1 = linspace(0, T, n1);
dt2 = 0.0043; n2 = size(numbound2, 1); t2 = linspace(0, T, n2);
dt3 = 0.0044; n3 = size(numbound3, 1); t3 = linspace(0, T, n3);

dt4 = 0.0026; n4 = size(numbound4, 1); t4 = linspace(0, T, n4);
dt5 = 0.0043; n5 = size(numbound5, 1); t5 = linspace(0, T, n5);
dt6 = 0.0044; n6 = size(numbound6, 1); t6 = linspace(0, T, n6);

figure();
set(gcf, 'Color', 'white');
subplot(1, 2, 1);
set(gca, 'FontSize', 16);
plot(t1, energy1, t2, energy2, t3, energy3);
title({'Mechanical energy of starcluster', '(\epsilon = 0.1)'});
xlabel('Time');
ylabel('Energy');
legend('N=100', 'N=250', 'N=1000');

subplot(1, 2, 2);
set(gca, 'FontSize', 16);
plot(t4, energy4, t5, energy5, t6, energy6);
title({'Mechanical energy of starcluster', '(\epsilon = 0.01)'});
xlabel('Time');
ylabel('Energy');
legend('N=100', 'N=250', 'N=1000');

figure();
set(gcf, 'Color', 'white');
subplot(1, 2, 1);
set(gca, 'FontSize', 16);
plot(t1, zeros(1, n1), 'k-', t1, virialenergy1, t2, virialenergy2, t3, virialenergy3);
title({'Virial energy energy in starcluster', '(\epsilon = 0.1)'});
xlabel('Time');
ylabel('Virial energy   2<T> + <V>');
legend('N=100', 'N=250', 'N=1000', 'Location', 'SE');

subplot(1, 2, 2);
set(gca, 'FontSize', 16);
plot(t4, zeros(1, n4), 'k-', t4, virialenergy4, t5, virialenergy5, t6, virialenergy6);
title({'Virial energy energy in starcluster', '(\epsilon = 0.01)'});
xlabel('Time');
ylabel('Virial energy   2<T> + <V>');
legend('N=100', 'N=250', 'N=1000', 'Location', 'SE');

figure();
set(gcf, 'Color', 'white');
subplot(1, 2, 1);
set(gca, 'FontSize', 16);
plot(t1, numbound1, t2, numbound2, t3, numbound3);
title({'Fraction of bound bodies vs. time', '(\epsilon = 0.1)'});
xlabel('Time');
ylabel('Fraction of bound bodies');
legend('N=100', 'N=250', 'N=1000', 'Location', 'SE');

subplot(1, 2, 2);
set(gca, 'FontSize', 16);
plot(t4, numbound4, t5, numbound5, t6, numbound6);
title({'Fraction of bound bodies vs. time', '(\epsilon = 0.01)'});
xlabel('Time');
ylabel('Fraction of bound bodies');
legend('N=100', 'N=250', 'N=1000');

%% Long T starcluster

T = 15;

numbound7 = load('../output/250_body_bound_15_1_0_0_0.10.dat')/250;
numbound8 = load('../output/250_body_bound_15_1_0_0_0.01.dat')/250;

n7 = size(numbound7, 1); dt7 = 0.0043; t7 = linspace(0, T, n7);
n8 = size(numbound8, 1); dt8 = 0.0043; t8 = linspace(0, T, n8);

figure();
set(gcf, 'Color', 'white');
set(gca, 'FontSize', 16);
plot(t7, numbound7, 'k', t8, numbound8);
title('Fraction of bound bodies vs. time for 250-body cluster');
xlabel('Time');
ylabel('Fraction of bound bodies');
legend('\epsilon = 0.1', '\epsilon = 0.01');


%% Part x - plotting radial distribution of particles. 

%positions = load('../output/1000_body_trajectories_6_1_0_0_0.10.dat');
%positions = positions(end, :);

N = size(positions, 2)/3;
r = zeros(1, N);
for i = 1:N
    r(i) = norm([positions(3*i-2) positions(3*i-1) positions(3*i)]);
end

numbins = 20;
[dist, centers] = hist(r, numbins);
half_width = (centers(2)-centers(1))/2;

density = zeros(1, numbins);

for i=1:numbins
    shell_volume = 4*pi*((centers(i)+half_width)^3 - (centers(i)-half_width)^3)/3;
    density(i) = dist(i)/shell_volume;
end

figure();
n_0 = 0.4;
r_0 = 7;
n = linspace(0, 100, 100);

n_r = n_0./(1+(n./r_0).^4);

plot(centers, density);
hold on;
plot(n, n_r, 'k');