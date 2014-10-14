
data = load('./trajectories/solarsystem12_trajectories.dat');
N = 12;

hold on;

for i = 1:N
    if i == 11
        plot(data(3*i-2,:), data(3*i-1,:), 'k');
    else
        plot(data(3*i-2,:), data(3*i-1,:));
    end
    
end

grid on;