% Animates the solar system or a part of the solar system
function [] = animate(positions, I, system, central, ax_size, im_size, dt) 

bodyimg = imread('../images/redbody.png');
n = size(positions, 1); % # of timesteps
set(gca, 'FontSize', 16);
set(gcf, 'Color', 'white');

for i = 5200
    plot(positions(1:i, I(1)), positions(1:i, I(2)), 'k');
    hold on;
    axis equal;
    for j = I
            plot(positions(1:i, 3*j-2), positions(1:i, 3*j-1));
    end
    
    for j = I
        image([positions(i, j*3-2)-im_size positions(i, j*3-2)+im_size],...
            [positions(i, j*3-1)+im_size positions(i, j*3-1)-im_size], bodyimg);
    end
    
    hold off;
    title({sprintf('The %s System', system), sprintf('%.1f weeks after 1.10.2014', dt*i)});   
    xlabel('Distance [AU]');
    ylabel('Distance [AU]');
    axis([positions(i, 3*central-2)-ax_size positions(i, 3*central-2)+ax_size...
        positions(i, 3*central-1)-ax_size positions(i, 3*central-1)+ax_size]); 
    pause(0.01);
    hold off;
end
end