%% Load data
data = load('./trajectories/solarsystem13_trajectories.dat');
N = 13;
n = size(data, 2);
AX = 5;
dt = 0.01;
bodyimg = imread('figures/redbody.png');

%% Static plot of entire system
hold on
for i = 1:N
    plot(data(6*i-5,:), data(6*i-4,:));
end
axis([-AX AX -AX AX]);

%% Animation of entire system
for i = 1:n
    plot(data(1, 1:i), data(2, 1:i), 'k');
    hold on;
    for j = 1:N
            plot(data(6*j-5,1:i), data(6*j-4,1:i));
    end
    %set(gca, 'Color', 'k');
    for j = 1:N
        h = 0.1;
        image([data(j*6-5, i)-h data(j*6-5, i)+h], [data(j*6-4, i)+h data(j*6-4, i)-h], bodyimg);
    end
    
    hold off;
    title(sprintf('The solar system %.1f weeks after 1.10.2014', dt*i));
    xlabel('Distance [AU]');
    ylabel('Distance [AU]');
    axis([-AX AX -AX AX]);
    pause(0.001);
    hold off;
end

%% Animation of Earth-Moon system
for i = 1:n
    plot(data(1, 1:i), data(2, 1:i), 'k');
    hold on;
    for j = [4 7]
        if j == 7
            plot(data(6*j-5,1:i), data(6*j-4,1:i), 'k');
        else
            plot(data(6*j-5,1:i), data(6*j-4,1:i));
        end
    end
    %set(gca, 'Color', 'k');
    
    for j = [4 7]
        h = 0.001;
        image([data(j*6-5, i)-h data(j*6-5, i)+h], [data(j*6-4, i)+h data(j*6-4, i)-h], bodyimg);
        %text(r(i,j*3-2) + h, r(i,j*3-1), planets{j}, 'Color', colors{j});
    end
    
    hold off;
    title(sprintf('The Earth-moon system %.1f weeks after 1.10.2014', dt*i));
    xlabel('Distance [AU]');
    ylabel('Distance [AU]');
    axis([data(6*4-5, i)-0.01 data(6*4-5, i)+0.01 data(6*4-4, i)-0.01 data(6*4-4, i)+0.01]); 
    pause(0.001);
    hold off;
end

%% Animation of Martian system
for i = 1:n
    a = 0.001;
    
    plot(data(1, 1:i), data(2, 1:i), 'k');
    hold on;
    for j = [5 8 9]
        if j == 8 || j == 9
            plot(data(6*j-5,1:i), data(6*j-4,1:i), 'k');
        else
            plot(data(6*j-5,1:i), data(6*j-4,1:i));
        end
    end
    %set(gca, 'Color', 'k');
    
    for j = [5 8 9]
        h = 0.00005;
        h2 = 0.00007;
        if j == 8 || j == 9
            image([data(j*6-5, i)-h data(j*6-5, i)+h], [data(j*6-4, i)+h data(j*6-4, i)-h], bodyimg);
        else
            image([data(j*6-5, i)-h2 data(j*6-5, i)+h2], [data(j*6-4, i)+h2 data(j*6-4, i)-h2], bodyimg);
        end
    end
    
    hold off;
    title(sprintf('The Mars system %.1f weeks after 1.10.2014', dt*i));
    xlabel('Distance [AU]');
    ylabel('Distance [AU]');
    axis([data(6*5-5, i)-a data(6*5-5, i)+a data(6*5-4, i)-a data(6*5-4, i)+a]); 
    pause(0.01);
    hold off;
end

%% Giovian system

for i = 1:n
    a = 0.05;
    
    plot(data(1, 1:i), data(2, 1:i), 'k');
    hold on;
    for j = [6  10 11 12 13]
        if j == 10 || j == 11 || j==12 ||j==13
            plot(data(6*j-5,1:i), data(6*j-4,1:i), 'k');
        else
            plot(data(6*j-5,1:i), data(6*j-4,1:i));
        end
    end
    %set(gca, 'Color', 'k');
    
    for j = [6 10 11 12 13]
        h = 0.0005;
        h2 = 0.0007;
        if j == 10 || j == 11 || j==12 ||j==13
            image([data(j*6-5, i)-h data(j*6-5, i)+h], [data(j*6-4, i)+h data(j*6-4, i)-h], bodyimg);
        else
            image([data(j*6-5, i)-h2 data(j*6-5, i)+h2], [data(j*6-4, i)+h2 data(j*6-4, i)-h2], bodyimg);
        end
    end
    
    hold off;
    title(sprintf('The Giovian system %.1f weeks after 1.10.2014', dt*i));
    xlabel('Distance [AU]');
    ylabel('Distance [AU]');
    axis([data(6*6-5, i)-a data(6*6-5, i)+a data(6*6-4, i)-a data(6*6-4, i)+a]); 
    pause(0.01);
  
    hold off;
end
