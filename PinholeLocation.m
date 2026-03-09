clc;
clear;
global all_clicked_points;
all_clicked_points = [];

data=load('E:\psf生成\vector_test\vector_test\PSF_10.mat');
psf=data.data;

[psfHeight, psfWidth, stack_depth] = size(psf);
rate=data.paraSim.zemit;



% r_min = 36200;
r_min = 0;
r_max = 44200;
spacing=250;
n = 32;
m =n * spacing / 2 / pi;

theta_min = 0;
theta_max = (r_max - r_min) / m;

theta = (theta_min : 0.001: theta_max);
r = r_min + m .* theta ;

coord = cell(n, 1);
d_theta = 2 * pi / n;
for idx = 1 : n
   coord{idx}(:,1) = r .* cos(theta - d_theta * (idx - 1));
   coord{idx}(:,2) = r .* sin(theta - d_theta * (idx - 1));
end
figure;
axis image tight;
axis xy;
hold on;
for idx = 1 : n
    plot(coord{idx}(:,1), coord{idx}(:,2));
end

arc_fun = @(theta) (1/(2*m))*((r_min + m.*theta).*sqrt((r_min + m.*theta).^2 + m.^2) ...
    + m^2.*log((r_min + m.*theta) + sqrt((r_min + m.*theta).^2 + m.^2))) ...
    - (1/(2*m))*(r_min*sqrt(r_min^2 + m^2) + m^2*log(r_min + sqrt(r_min^2 + m^2)));
total_length = arc_fun(theta_max);
target_lengths = 0:spacing:total_length;
theta_samples = zeros(size(target_lengths));

for i = 1:length(target_lengths)
    target = target_lengths(i);
    theta_samples(i) = fzero(@(th) arc_fun(th) - target, [0, theta_max]);
end
r_samples = r_min + m .* theta_samples;

pinhole_coords = cell(n,1);
for idx = 1:n
    
    pinhole_coords{idx}(:,1)=r_samples.*cos(theta_samples - d_theta * (idx - 1));
    pinhole_coords{idx}(:,2)=r_samples.*sin(theta_samples - d_theta * (idx - 1));
       

end

fig=figure;
hold on;
h = gobjects(n,1);
for idx = 1 : n
    h(idx)=scatter(pinhole_coords{idx}(:,1), pinhole_coords{idx}(:,2));
end

set(fig, 'WindowButtonDownFcn', @recordClick);
for i = 1:numel(h)
    set(h(i), 'ButtonDownFcn', @recordClick);
end

axis image tight;
axis xy;

drawnow; 
pause(0.1);

    function recordClick(~, evt)
        
        global all_clicked_points;      

        pos = evt.IntersectionPoint(1:2);
        all_clicked_points = [all_clicked_points; pos];
        
        scatter(pos(1), pos(2), 30, 'r', 'filled');
        
        num_points = size(all_clicked_points, 1);
        fprintf('Point %d saved: (%.2f, %.2f)\n', num_points, pos(1), pos(2));
        
        if num_points == 2
            save('clicked_points.mat', 'all_clicked_points');
        elseif num_points > 2
             save('clicked_points.mat', 'all_clicked_points');
        end
    end
