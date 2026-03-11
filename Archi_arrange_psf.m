clc;
clear;

position=load("pinhole_coords.mat");
coords=position.all_coords;

data=load('E:\psf生成\vector_test\vector_test\PSF_10.mat');
psf=data.data;
rate=data.paraSim.zemit;

clicked_points=load("clicked_points.mat");
center_pinhole_coord=clicked_points.all_clicked_points;
center_pinhole_coord(1,:)=[];
[center_pinhole_coords, ~] = unique(center_pinhole_coord, 'rows', 'stable');
pixel_size=50e-9;
pinhole_d=0.5e-6;

[psfHeight, psfWidth, stack_depth] = size(psf);
circle_center = zeros(2*psfHeight, 2*psfWidth);
center_x = round(psfWidth );
center_y = round(psfHeight );
r=pinhole_d/(pixel_size*2);


for y = 1:2*psfHeight
    for x = 1:2*psfWidth
        distance = sqrt((x - center_x)^2 + (y - center_y)^2);
        if distance <= r
            circle_center(y, x) = 1;
        end
    end
end

crosstalk_num=size(center_pinhole_coords, 1);
figure('Color', 'w');
hold on;
grid on;
colors = lines(crosstalk_num);
legend_str = cell(crosstalk_num, 1);
for k=1:crosstalk_num
    crosstalk_pinhole=[];
    dis=[];
    current_center=center_pinhole_coords(k,:);
    for i = 1:size(coords, 1)
        dist = sqrt((coords(i,1) - current_center(1))^2 + ...
                    (coords(i,2) - current_center(2))^2);
        
        if dist <= 200*r&&dist>1
            dis=[dis,dist];
            crosstalk_pinhole = [crosstalk_pinhole; coords(i,:)];  
        end
    end
    
    % mask=zeros(2*psfHeight,2*psfWidth);
    pinhole_circle_center = zeros(length(crosstalk_pinhole),2*psfHeight, 2*psfWidth);
    crosstalk_pinhole_center_x=zeros(length(crosstalk_pinhole),1);
    crosstalk_pinhole_center_y=zeros(length(crosstalk_pinhole),1);
    mask=circle_center;
    for j=1:length(crosstalk_pinhole)
        crosstalk_pinhole_center_x(j)=(crosstalk_pinhole(j,1)-current_center(1,1))/1e8/pixel_size+center_x;
        crosstalk_pinhole_center_y(j)=(crosstalk_pinhole(j,2)-current_center(1,2))/1e8/pixel_size+center_y;
        for y = 1:2*psfHeight
            for x = 1:2*psfWidth
                distance = sqrt((x - crosstalk_pinhole_center_x(j))^2 + (y - crosstalk_pinhole_center_y(j))^2);
                if distance <= r
                    pinhole_circle_center(j,y,x) = 1;
                end
            end
        end
        mask=mask+squeeze(pinhole_circle_center(j, :, :));
    end
    mask_start_conv=conv2(mask,mask,"same");
    mask_sum=mask_start_conv;
    angles = 1:1:359;
    for a = angles
        mask_rot = imrotate(mask, a, 'bilinear', 'crop');
        mask_conv=conv2(mask_rot,mask_rot,"same");
        mask_sum=mask_sum+mask_conv;
    end
    mask_avg=mask_sum/(length(angles)+1);
    cutoff_pix = (pinhole_d * 20) / (pixel_size);
    [X, Y] = meshgrid(1:2*psfWidth, 1:2*psfHeight);
    mask_avg(sqrt((X-center_x).^2 + (Y-center_y).^2) > cutoff_pix) = 0;
    for depth=1:stack_depth
        psf_exc(:,:,depth)= conv2(psf(:,:,depth) , mask_avg,'same');
    end
    psf_overall=psf.*psf_exc;
    for depth=1:stack_depth
        I_psf(depth)=sum(psf_overall(:,:,depth),"all");
    end
    I_psf=I_psf./max(I_psf);
    plot(1:stack_depth, I_psf, 'LineWidth', 1.5, 'Color', colors(k,:));
    legend_str{k} = sprintf('Center Point %d', k);
    
    drawnow;
end
d=pinhole_d*5/pixel_size;
a = d;             
pinhole_num=200*r/a;        
maxLayer = pinhole_num;     

e1 = [1, 0];
e2 = [0.5, sqrt(3)/2];

local_points = [];
for n = 0:maxLayer
    
    if n == 0
        local_points = [local_points; 0 0];   
        continue;
    end
    
    x = n;
    y = 0;
    pos = x*e1 + y*e2;
    
    dirs = [ 0 -1;   -1 0;   -1 1;  0 1;   1 0;   1 -1 ];
    
    for side = 1:6
        for step = 1:n
            local_points = [local_points; pos];
            x = x + dirs(side,1);
            y = y + dirs(side,2);
            pos = x*e1 + y*e2;
        end
    end
end
points = local_points * a ;
points = unique(points, 'rows');
simu_center=zeros(length(points),2*psfHeight,2*psfWidth);
mask_simu=circle_center;
for i=1:size(points,1)
    simu_x = center_x + points(i,1);
    simu_y = center_y + points(i,2);
    for y = 1:2*psfHeight
        for x = 1:2*psfWidth
            distance = sqrt((x - simu_x)^2 + (y - simu_y)^2);
            if distance <= r
                simu_center(i,y,x) = 1;
            end
        end
    end
    mask_simu=mask_simu+squeeze(simu_center(i, :, :));
end
mask_simu_conv=conv2(mask_simu,mask_simu,"same");
mask_simu_sum=mask_simu_conv;
for a = angles
    mask_simu_rot = imrotate(mask_simu, a, 'bilinear', 'crop');
    mask_simu_conv=conv2(mask_simu_rot,mask_simu_rot,"same");
    mask_simu_sum=mask_simu_sum+mask_simu_conv;
end
mask_simu_avg=mask_simu_sum/(length(angles)+1);
mask_simu_avg(sqrt((X-center_x).^2 + (Y-center_y).^2) > cutoff_pix) = 0;
for depth=1:stack_depth
    psf_simu_exc(:,:,depth)= conv2(psf(:,:,depth) , mask_avg,'same');
end
psf_simu_overall=psf.*psf_simu_exc;
for depth=1:stack_depth
    I_simu_psf(depth)=sum(psf_simu_overall(:,:,depth),"all");
end
I_simu_psf=I_simu_psf./max(I_simu_psf);
plot(1:stack_depth, I_simu_psf, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
figure;
imagesc(mask_simu_avg);
colormap jet;
axis image tight;
axis xy;