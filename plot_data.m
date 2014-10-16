%% Load data
data = load('./trajectories/solarsystem18_trajectories.dat');
N = 18;
n = size(data, 2);
AX = 30;
dt = 0.05;
bodyimg = imread('redbody.png');

%% calcs
hold on
for i = 1:N
    if i == 11
        plot(data(3*i-2,:), data(3*i-1,:), 'k');
    else
        plot(data(3*i-2,:), data(3*i-1,:));
    end
    axis([-AX AX -AX AX]);
end

%% plot 3d

hold on
for i = 1:N
   plot3(data(3*i-2,:), data(3*i-1,:), data(3*i,:));
end


%% moon earth animation
AX = 0.01;
for i = 1:n
    
    plot(data(1, 1:i), data(2, 1:i), 'k');
    hold on;
    for j = [4 11]
        if j == 11
            plot(data(3*j-2,1:i), data(3*j-1,1:i), 'k');
        else
            plot(data(3*j-2,1:i), data(3*j-1,1:i));
        end
    end
    %set(gca, 'Color', 'k');
    
    for j = [4 11]
        h = 0.001;
        image([data(j*3-2, i)-h data(j*3-2, i)+h], [data(j*3-1, i)+h data(j*3-1, i)-h], bodyimg);
        %text(r(i,j*3-2) + h, r(i,j*3-1), planets{j}, 'Color', colors{j});
    end
    
    hold off;
    title(sprintf('The Earth-moon system %.1f weeks after 1.10.2014', dt*i));
    xlabel('Distance [AU]');
    ylabel('Distance [AU]');
    axis([data(3*4-2, i)-0.05 data(3*4-2, i)+0.05 data(3*4-1, i)-0.05 data(3*4-1, i)+0.05]); 
    %%legend(planets(1:no_planets), 'TextColor', 'white');
    pause(0.001);
  
    %%film = addframe(film, getframe(fig));
    hold off;
end

%% Mars system animation
for i = 1:n
    a = 0.001;
    
    plot(data(1, 1:i), data(2, 1:i), 'k');
    hold on;
    for j = [5 12 13]
        if j == 12 || j == 13
            plot(data(3*j-2,1:i), data(3*j-1,1:i), 'k');
        else
            plot(data(3*j-2,1:i), data(3*j-1,1:i));
        end
    end
    %set(gca, 'Color', 'k');
    
    for j = [5 12 13]
        h = 0.00005;
        h2 = 0.00007;
        if j == 12 || j == 13
            image([data(j*3-2, i)-h data(j*3-2, i)+h], [data(j*3-1, i)+h data(j*3-1, i)-h], bodyimg);
        else
            image([data(j*3-2, i)-h2 data(j*3-2, i)+h2], [data(j*3-1, i)+h2 data(j*3-1, i)-h2], bodyimg);
        end
    end
    
    hold off;
    title(sprintf('The Mars system %.1f weeks after 1.10.2014', dt*i));
    xlabel('Distance [AU]');
    ylabel('Distance [AU]');
    axis([data(3*5-2, i)-a data(3*5-2, i)+a data(3*5-1, i)-a data(3*5-1, i)+a]); 
    %%legend(planets(1:no_planets), 'TextColor', 'white');
    pause(0.01);
  
    %%film = addframe(film, getframe(fig));
    hold off;
end

%% Giovian system animation
for i = 1:n
    a = 0.05;
    
    plot(data(1, 1:i), data(2, 1:i), 'k');
    hold on;
    for j = [6  14 15 16 17]
        if j == 14 || j == 15 || j==16 ||j==17
            plot(data(3*j-2,1:i), data(3*j-1,1:i), 'k');
        else
            plot(data(3*j-2,1:i), data(3*j-1,1:i));
        end
    end
    %set(gca, 'Color', 'k');
    
    for j = [6 14 15 16 17]
        h = 0.0005;
        h2 = 0.0007;
        if j == 14 || j == 15 || j==16 ||j==17
            image([data(j*3-2, i)-h data(j*3-2, i)+h], [data(j*3-1, i)+h data(j*3-1, i)-h], bodyimg);
        else
            image([data(j*3-2, i)-h2 data(j*3-2, i)+h2], [data(j*3-1, i)+h2 data(j*3-1, i)-h2], bodyimg);
        end
    end
    
    hold off;
    title(sprintf('The Mars system %.1f weeks after 1.10.2014', dt*i));
    xlabel('Distance [AU]');
    ylabel('Distance [AU]');
    axis([data(3*6-2, i)-a data(3*6-2, i)+a data(3*6-1, i)-a data(3*6-1, i)+a]); 
    %%legend(planets(1:no_planets), 'TextColor', 'white');
    pause(0.01);
  
    %%film = addframe(film, getframe(fig));
    hold off;
end

%% Entire solarsystem animation
for i = 1:n
    
    plot(data(1, 1:i), data(2, 1:i), 'k');
    hold on;
    for j = 1:N
            plot(data(3*j-2,1:i), data(3*j-1,1:i));
    end
    %set(gca, 'Color', 'k');
    
%     for j = 1:no_planets
%         h = im_sizes(j);
%         image([r(i,j*3-2)-h r(i,j*3-2)+h], [r(i,j*3-1)+h r(i,j*3-1)-h], images{j});
%         %text(r(i,j*3-2) + h, r(i,j*3-1), planets{j}, 'Color', colors{j});
%     end
    
    hold off;
    title(sprintf('The solar system %.1f weeks after 1.10.2014', dt*i));
    xlabel('Distance [AU]');
    ylabel('Distance [AU]');
    axis([-AX, AX, -AX, AX]);
    %%legend(planets(1:no_planets), 'TextColor', 'white');
    pause(0.001);
  
    %%film = addframe(film, getframe(fig));
    hold off;
end

grid on;