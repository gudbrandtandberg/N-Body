%% first attempt at plotting open cluster

data = load('./trajectories/cluster100_trajectories.dat');
N = 100;
n = size(data, 2);

for i = 1:n
    hold on;
    
    for j = 1:N
        plot3(data(6*j-5,1:i), data(6*j-4, 1:i), data(6*j-3, 1:i));
        
    end
    
    hold off;
    pause(0.001);
end

