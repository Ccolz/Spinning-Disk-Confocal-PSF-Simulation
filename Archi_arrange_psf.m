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
colors = lines(crosstalk_num); % 生成不同颜色的线条
legend_str = cell(crosstalk_num, 1);
for k=1:crosstalk_num
    crosstalk_pinhole=[];
    dis=[];
    current_center=center_pinhole_coords(k,:);
    for i = 1:size(coords, 1)
        dist = sqrt((coords(i,1) - current_center(1))^2 + ...
                    (coords(i,2) - current_center(2))^2);
        
        if dist <= pinhole_d*2*1e9&&dist>1
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
    angles = 1:1:360;
    for a = angles
        mask_rot = imrotate(mask, a, 'bilinear', 'crop');
        mask_conv=conv2(mask_rot,mask_rot,"same");
        mask_sum=mask_sum+mask_conv;
    end
    mask_avg=mask_sum/length(angles);
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


% figure;
% imagesc(mask_avg);
% colormap jet;
% axis image tight;
% axis xy;