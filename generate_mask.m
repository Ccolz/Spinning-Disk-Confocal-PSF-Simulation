clc;
clear;

data=load('.\PSF_10.mat');
psf=data.data;

[psfHeight, psfWidth, stack_depth] = size(psf);

pixel_size=50e-9;
tmp_z=61;
z_position=2000-50*(tmp_z-1);

t_list=[1,2,3,5,6];
% t_list=[1];
pinhole_d_list=[0.5e-6];

circle_center = zeros(4*psfHeight, 4*psfWidth);
sd_circle_center = zeros(4*psfHeight, 4*psfWidth);

centers=zeros(4*psfHeight, 4*psfWidth);
center_x = round(psfWidth*2 );
center_y = round(psfHeight*2 );
z_axis=zeros(1,stack_depth);
I_WF_z=zeros(1,stack_depth);
I_con_z=zeros(length(pinhole_d_list),stack_depth);
I_sd_z=zeros(length(t_list),length(pinhole_d_list),stack_depth);
mask_results = cell(length(pinhole_d_list), length(t_list));

for i_d=1:length(pinhole_d_list)
    for i_t=1:length(t_list)
        pinhole_d=pinhole_d_list(i_d);
        t=t_list(i_t);
        r=pinhole_d/(pixel_size*2);
        d=2*r*t;
        %%Pinhole of Confocal
        for y = 1:4*psfHeight
            for x = 1:4*psfWidth
                distance = sqrt((x - center_x)^2 + (y - center_y)^2);
                if distance <= r
                    circle_center(y, x) = 1;
                end
            end
        end
        pinhole_intensity=conv2(circle_center,circle_center,'same');
        border=10*pinhole_d/pixel_size;
        pinhole_num=floor(border/d);
        %%Mask of Spinning Disk
        %%Arrangement of the points in the real mask
          
        a = d;             
        
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
        points = local_points * a + [center_x, center_y];
        points = unique(points, 'rows');
        mask0=zeros(4*psfHeight, 4*psfWidth);
        mask=zeros(4*psfHeight, 4*psfWidth);

        for i=1:length(points)
            sd_circle_center(:) = 0;
            
             cx = round(points(i,1));
             cy = round(points(i,2));

    
            mask(cy, cx) = mask(cy, cx) + 1; %the original points 
        end
        %%rotate the mask
        angles = 0:1:60; 
        mask_sum = zeros(4*psfHeight, 4*psfWidth);
        for a = angles 
            mask_rot = imrotate(mask, a, 'bilinear', 'crop');

            mask_conv=conv2(mask_rot,mask_rot,"same"); %convolution of points with different angles
            
            mask_bin = mask_conv ~= 0;
            [cy_list, cx_list] = find(mask_bin);
            %%rearrange the convolution points
            conv_points=[cx_list,cy_list];
            merged = false(size(conv_points,1),1);

            new_points = [];
            cx_new = [];
            cy_new = [];
            d_merge = r/2;
            for i = 1:size(conv_points,1)

                if merged(i)
                    continue;
                end
            
                dist = sqrt( (conv_points(:,1)-conv_points(i,1)).^2 + ...
                             (conv_points(:,2)-conv_points(i,2)).^2 );
            
                idx = dist < d_merge;
            
                cx = round(mean(conv_points(idx,1)));
                cy = round(mean(conv_points(idx,2)));
            
                cx_new(end+1) = cx;
                cy_new(end+1) = cy;
                new_points(end+1,:) = [cx, cy];
            
                merged(idx) = true;
            end

            %%replace the points with the result of single-point
            %%convolution
            mask_sd = zeros(size(mask_bin));

            [H, W] = size(pinhole_intensity);
            cx0 = round(W/2);
            cy0 = round(H/2);

            for k = 1:length(cx_new)

                cx = cx_new(k);
                cy = cy_new(k);

                x1 = cx - cx0 + 1;
                x2 = cx - cx0 + W;
                y1 = cy - cy0 + 1;
                y2 = cy - cy0 + H;

                xs = max(1, x1);  xe = min(size(mask_sd,2), x2);
                ys = max(1, y1);  ye = min(size(mask_sd,1), y2);

                px1 = xs - x1 + 1;
                px2 = px1 + (xe - xs);
                py1 = ys - y1 + 1;
                py2 = py1 + (ye - ys);

                mask_sd(ys:ye, xs:xe) = mask_sd(ys:ye, xs:xe) + ...
                                        pinhole_intensity(py1:py2, px1:px2);
            end

            mask_sum = mask_sum + mask_sd;
            
        end 
            
        mask_avg = mask_sum / length(angles);
        % mask_results = cell(length(pinhole_d_list), length(t_list)); 
        % mask_results{i_d, i_t} = mask_avg;
        %% 
        disk_psf_exc_avg=zeros(psfHeight,psfWidth,stack_depth);
        for depth_Idx=1:stack_depth
            disk_psf_exc_avg(:,:,depth_Idx)= conv2(psf(:,:,depth_Idx) , mask_avg,'same'); %Spinning Disk Confocal Excitation PSF
        end
        disk_psf_avg=disk_psf_exc_avg.*psf; %Spinning Disk Confocal PSF
        con_psf_det=zeros(psfHeight,psfWidth,stack_depth);
        for depth_Idx=1:stack_depth
            con_psf_det(:,:,depth_Idx)= conv2(psf(:,:,depth_Idx) , circle_center,'same');%Confocal Detection PSF
        end
        con_psf=con_psf_det.*psf; %Confocal PSF

        wf_psf=psf; %Wide Field PSF

        for depth=1:stack_depth
            I_sd_z(i_t,i_d,depth)=sum(disk_psf_avg(:,:,depth),"all");
        end
        I_sd_z(i_t,i_d,:)=(I_sd_z(i_t,i_d,:))/(max(I_sd_z(i_t,i_d,:)));
    end
    for depth=1:stack_depth
        I_con_z(i_d,depth)=sum(con_psf(:,:,depth),"all");
    end
    I_con_z(i_d,:)=(I_con_z(i_d,:))/(max(I_con_z(i_d,:)));
end
for depth=1:stack_depth
    I_WF_z(depth)=sum(wf_psf(:,:,depth),"all");
    z_axis(depth)=2000-50*(depth-1);
end
I_WF_z=(I_WF_z)/(max(I_WF_z));

figure;
colors = lines(length(t_list));
for i_d = 1:length(pinhole_d_list)
    subplot(1, length(pinhole_d_list), i_d);
    axis tight;
    hold on;
    plot(z_axis,squeeze(I_WF_z), 'k--', 'LineWidth', 1.5, 'DisplayName', 'WF');
    plot(z_axis,squeeze(I_con_z(i_d,:)), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Confocal');
    for i_t = 1:length(t_list)
        sd = squeeze(I_sd_z(i_t, i_d, :));
        plot(z_axis,sd, ':', 'Color', colors(i_t,:), 'LineWidth', 1.5, 'DisplayName', sprintf('SD distance=%dd', t_list(i_t)));
    end

    hold off;
    title(sprintf('Pinhole diameter = %.2f μm', pinhole_d_list(i_d)*1e8), 'FontSize', 12);
    xlabel('z (nm)');
    ylabel('Normalized intensity');
    legend('show', 'Location', 'northeastoutside');
    grid on;
end
figure;
imagesc(mask);
title(sprintf('Pinhole diameter = %.2f μm, distance = %.2f μm', pinhole_d_list(i_d)*1e8, t_list(i_t)*pinhole_d_list(i_d)*1e8), 'FontSize', 12);

colormap jet;
axis image tight;
axis xy;
save('sd_psf.mat', 'disk_psf_avg');
save('con_psf.mat','con_psf')
