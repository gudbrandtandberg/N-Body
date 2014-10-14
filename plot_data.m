
data = load('./solarsystem_trajectories.dat');
N = 9;

hold on;

for i = 1:N
    plot(data(3*i-2,:), data(3*i-1,:));
end

grid on;