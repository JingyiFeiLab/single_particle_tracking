%clear all
close all

sample = 'C4';
saveFile = strcat(['/Users/reyer/Documents/MATLAB/SOURCE_CODES/sample_images_matt/Matt_Microscope/Pol2_31521/',sample,'/']);

repositionFigures = 1; % 1 to reposition/resize figures according to my laptop, 0 to leave as is
saveFigures = 1; % 1 to save figures to "saveFile" location, 0 to not save (automatically)
low_pass_register = 0; % 1 to register the low pass images, 0 to register original images


Npts = 2;
Eps = 15;

gfp_path=strcat(['/Users/reyer/Documents/MATLAB/SOURCE_CODES/sample_images_matt/Matt_Microscope/Pol2_31521/',sample,'/',sample,'_M3_488.tif']);
gfp_post_path=strcat(['/Users/reyer/Documents/MATLAB/SOURCE_CODES/sample_images_matt/Matt_Microscope/Pol2_31521/',sample,'/',sample,'_M3_647_Pol2.tif']);

fixed_im_path=strcat(['/Users/reyer/Documents/MATLAB/SOURCE_CODES/sample_images_matt/Matt_Microscope/Pol2_31521/',sample,'/',sample,'_STORM.png']);
moving_im_path=strcat(['/Users/reyer/Documents/MATLAB/SOURCE_CODES/sample_images_matt/Matt_Microscope/Pol2_31521/',sample,'/',sample,'_sptPALM.png']);

storm_file = strcat(['/Users/reyer/Documents/MATLAB/SOURCE_CODES/sample_images_matt/Matt_Microscope/Pol2_31521/',sample,'/',sample,'_STORM.csv']);
msd_file = (['/Users/reyer/Documents/MATLAB/SOURCE_CODES/sample_images_matt/Matt_Microscope/Pol2_31521/',sample,'/onestep_MSD_list_min_3.dat']);

storm = csvread(storm_file,1,0);
msd = fopen(msd_file);
msd = textscan(msd,'%s %s %s %s %s');
msd = [str2double(msd{1}) str2double(msd{2}) str2double(msd{3}) str2double(msd{4}) str2double(msd{5})];

storm_pix_size = 173; %Microns
msd_pix_size = .130;

int_thresh = .0001; % Intensity Threshold
shape2D_thresh = 2;  % <--- Splitting threshold, you can change this % higher is more generous   (default=5.5)  8

background_thresh = .2; % lower is generous (default=.1~.2) .1
pix_neigh = 5;%floor((.08/pix_size)*12);% <- If you have no idea, try this. Default=8 or 11  (5 worked well) 25

low_pass_thresh = 0.025; % less is more generous (try 0.01 or 0.02 )
high_pass_thresh = 0.075;

cluster_file = sprintf('%s_Cluster_Centers_ColorCode%d_%d_%d.dat',output,ColorCode,Npts,Eps);

total_centers = TotalCenters_master(cluster_file,Start_Npts,End_Npts,Increment_Npts,Start_Eps,End_Eps,Increment_Eps);

storm_x_col = 2; % Column for x coordinates in STORM csv
storm_y_col = 3; % Column for y coordinates in STORM csv

cluster_x_col = 1; % Column for x coordinates in total_centers
cluster_y_col = 2; % Column for y coordinates in total_centers
cluster_radius_col = 4; % Column for cluster radius in total_centers

msd_x_col = 1; % Column for x coordinates in MSD .dat
msd_y_col = 2; % Column for x coordinates in MSD .dat

% Path to main file (i.e. channel) that you will use for segmentation
%filepath_dic = strcat(['/Users/reyer/Documents/MATLAB/SOURCE_CODES/sample_images_matt/Matt_Microscope/August_24_17_convert/+SgrS/t20/dic',num2str(cell_num),'.tif']);
filepath_dic = gfp_path;
stack_o = mat2gray(imread(filepath_dic));
stack_2 = mat2gray(imread(gfp_post_path));

stack_fixed = mat2gray(imread(gfp_path));
stack_moving = mat2gray(imread(gfp_post_path));

X = size(stack_fixed,1) ;                              % CapitalX is the number of X pixels in the reference image
Y = size(stack_fixed,2) ; 

se = [1 1 1; 1 1 1 ; 1 1 1]; % Structuring Element for basic Erosion and dilation

fixed_low_pass = low_pass(stack_fixed,low_pass_thresh); %.025
fixed_low_pass = fixed_low_pass./max(fixed_low_pass(:));
moving_low_pass = low_pass(stack_moving,low_pass_thresh); %.025
moving_low_pass = moving_low_pass./max(moving_low_pass(:));

[optimizer,metric] = imregconfig('multimodal');
if low_pass_register == 1
    stack_recovered = imregister(moving_low_pass,fixed_low_pass,'rigid',optimizer,metric);
    tform = imregtform(moving_low_pass,fixed_low_pass,'rigid', optimizer, metric);
    register_low_pass = imwarp(moving_low_pass,tform,'OutputView',imref2d(size(fixed_low_pass)));
    register_storm = imwarp(stack_moving,tform,'OutputView',imref2d(size(fixed_low_pass)));
else
    stack_recovered = imregister(stack_moving,stack_fixed,'rigid',optimizer,metric);
    tform = imregtform(stack_moving,stack_fixed,'rigid', optimizer, metric);
    register_low_pass = imwarp(stack_moving,tform,'OutputView',imref2d(size(fixed_low_pass)));
    register_storm = imwarp(stack_moving,tform,'OutputView',imref2d(size(fixed_low_pass)));
end

C = imfuse(fixed_low_pass,moving_low_pass);%figure(3);imshow(C) Pre-registration overlap
C2 = imfuse(fixed_low_pass,stack_recovered);%figure(4);imshow(C2) Post-registration overlap

% [~,fixed_threshold] = edge(stack_fixed,'sobel');
% fudgeFactor = 0.5;
% BW_fixed = edge(stack_fixed,'sobel',fixed_threshold * fudgeFactor);
% 
% [~,moving_threshold] = edge(stack_moving,'sobel');
% fudgeFactor = 0.5;
% BW_moving = edge(stack_moving,'sobel',moving_threshold * fudgeFactor);



storm_x = storm(:,storm_x_col)./storm_pix_size;
storm_y = storm(:,storm_y_col)./storm_pix_size;

msd_x = msd(:,msd_x_col)./msd_pix_size;
msd_y = msd(:,msd_y_col)./msd_pix_size;
msd_coef = sqrt(msd(:,4)).*1000;

cluster_x = total_centers(:,cluster_x_col)./storm_pix_size;
cluster_y = total_centers(:,cluster_y_col)./storm_pix_size;
cluster_rad = total_centers(:,cluster_radius_col)./storm_pix_size;

figure(1);
imshow(5*stack_fixed);hold on;scatter(msd_y,msd_x,5,'filled','r')
continuetranslation=1;
while continuetranslation==1
    prompt={'X_Translation(- = Down, + = Up):','Y_Translation(- = right, + = left):','Continue the translation(1 == Yes && 0== NO) ?'};   % A box will take in the values for the X/Ytranslation
    mask_title='Translation';                             % The title of the box
    answer=inputdlg(prompt,mask_title);
    Xtranslation = str2num(answer{1}); 
    Ytranslation = str2num(answer{2});
    if isempty(answer{1})
        Xtranslation = 0;
    end

    if isempty(answer{2})
        Ytranslation = 0;
    end
    continuetranslation = str2num(answer{3});
    msd_x=msd_x-Xtranslation;                                % Translated
    msd_y=msd_y-Ytranslation;
    imshow(stack_fixed);hold on;scatter(msd_y,msd_x,5,'filled','r')
end

title('MSD Spots Over palm 488')
legend('PALM')
xlabel('Cell X Axis')
ylabel('Cell Y Axis')
if repositionFigures == 1
    set(gcf,'position',[2,901,386,319])
elseif repositionFigures == 2
    set(gcf,'position',[6,486,386,319])
end

if saveFigures == 1
    set(gcf,'PaperPositionMode','auto')
    file1 = strcat([saveFile,'/msd_over_488']);
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([saveFile,'/msd_over_488.fig']);
    savefig(gcf,file1_fig)
end

figure(2);
imshow(5*stack_fixed);hold on;scatter(storm_x,storm_y,1,'filled','r');hold on;scatter(msd_y,msd_x,5,'filled','g');
title('Spots Overlap, Original')
legend('STORM','PALM')
xlabel('Cell X Axis')
ylabel('Cell Y Axis')
if repositionFigures == 1
    set(gcf,'position',[116,1106,676,597])
elseif repositionFigures == 2
    set(gcf,'position',[393,207,676,597])
end

if saveFigures == 1
    set(gcf,'PaperPositionMode','auto')
    file1 = strcat([saveFile,'/spots_overlap']);
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([saveFile,'/spots_overlap.fig']);
    savefig(gcf,file1_fig)
end

[storm_x_reg,storm_y_reg] = transformPointsForward(tform,storm_x,storm_y);
[cluster_x_reg,cluster_y_reg] = transformPointsForward(tform,cluster_x,cluster_y);

figure(3)
imshow(5*stack_fixed);hold on;scatter(storm_x_reg,storm_y_reg,1,'filled','r');hold on;scatter(msd_y,msd_x,5,'filled','g');
title('Spots Overlap, Registered')
legend('STORM','PALM')
xlabel('Cell X Axis')
ylabel('Cell Y Axis')
if repositionFigures == 1
    set(gcf,'position',[796,1106,676,597])
elseif repositionFigures == 2
    set(gcf,'position',[765,208,676,597])
end
if saveFigures == 1
    set(gcf,'PaperPositionMode','auto')
    file1 = strcat([saveFile,'/spots_overlap_register']);
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([saveFile,'/spots_overlap_register.fig']);
    savefig(gcf,file1_fig)
end

figure(4);scatter(storm_x,storm_y,1,'filled','r');hold on;scatter(msd_y,msd_x,5,'filled','g');
title('Spots Overlap, Original')
legend('STORM','PALM')
xlabel('Cell X Axis')
ylabel('Cell Y Axis')
if repositionFigures == 1
    set(gcf,'position',[367,904,556,421])
elseif repositionFigures == 2
    set(gcf,'position',[1,1,556,421])
end
if saveFigures == 1
    set(gcf,'PaperPositionMode','auto')
    file1 = strcat([saveFile,'/spots_overlap_scatter']);
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([saveFile,'/spots_overlap_scatter.fig']);
    savefig(gcf,file1_fig)
end

figure(5);scatter(storm_x_reg,storm_y_reg,1,'filled','r');hold on;scatter(msd_y,msd_x,5,'filled','g');
title('Spots Overlap, Original')
legend('STORM','PALM')
xlabel('Cell X Axis')
ylabel('Cell Y Axis')
if repositionFigures == 1
    set(gcf,'position',[924,904,556,421])
elseif repositionFigures == 2
    set(gcf,'position',[558,1,556,421])
end
if saveFigures == 1
    set(gcf,'PaperPositionMode','auto')
    file1 = strcat([saveFile,'/spots_overlap_scatter']);
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([saveFile,'/spots_overlap_scatter.fig']);
    savefig(gcf,file1_fig)
end


storm_x_coord = int32(storm_x_reg);
storm_y_coord = int32(storm_y_reg);

storm_x_coord_nor = int32(storm_x);
storm_y_coord_nor = int32(storm_y);

storm_x_coord(storm_x_coord<1) = 1;
storm_y_coord(storm_y_coord<1) = 1;
storm_x_coord(storm_x_coord>=X) = X;
storm_y_coord(storm_y_coord>=Y) = Y;

storm_x_coord_nor(storm_x_coord_nor<1) = 1;
storm_y_coord_nor(storm_y_coord_nor<1) = 1;
storm_x_coord_nor(storm_x_coord_nor>=X) = X;
storm_y_coord_nor(storm_y_coord_nor>=Y) = Y;

msd_x_coord = int32(msd_x);
msd_y_coord = int32(msd_y);

msd_x_coord(msd_x_coord<1) = 1;
msd_y_coord(msd_y_coord<1) = 1;
msd_x_coord(msd_x_coord>=X) = X;
msd_y_coord(msd_y_coord>=Y) = Y;

cluster_x_coord = int32(cluster_x_reg);
cluster_y_coord = int32(cluster_y_reg);

cluster_x_coord_nor = int32(cluster_x);
cluster_y_coord_nor = int32(cluster_y);

cluster_x_coord(cluster_x_coord<1) = 1;
cluster_y_coord(cluster_y_coord<1) = 1;
cluster_x_coord(cluster_x_coord>=X) = X;
cluster_y_coord(cluster_y_coord>=Y) = Y;

cluster_x_coord_nor(cluster_x_coord_nor<1) = 1;
cluster_y_coord_nor(cluster_y_coord_nor<1) = 1;
cluster_x_coord_nor(cluster_x_coord_nor>=X) = X;
cluster_y_coord_nor(cluster_y_coord_nor>=Y) = Y;

msd_im = zeros(size(stack_fixed));
diffusion_im = zeros(size(stack_fixed));
storm_im = zeros(size(stack_fixed));
storm_im_nor = zeros(size(stack_fixed));

cluster_im = zeros(size(stack_fixed));
cluster_im_nor = zeros(size(stack_fixed));

for i = 1:length(storm_x_reg)
        
        storm_im(storm_y_coord(i),storm_x_coord(i)) = storm_im(storm_y_coord(i),storm_x_coord(i)) + 1;
        storm_im_nor(storm_y_coord_nor(i),storm_x_coord_nor(i)) = storm_im_nor(storm_y_coord_nor(i),storm_x_coord_nor(i)) + 1;
        
end

for i = 1:length(msd)
        
        msd_im(msd_x_coord(i),msd_y_coord(i)) = msd_im(msd_x_coord(i),msd_y_coord(i)) + 1;
        diffusion_im(msd_x_coord(i),msd_y_coord(i)) = diffusion_im(msd_x_coord(i),msd_y_coord(i)) + msd_coef(i);
end

for i = 1:length(cluster_x_reg)
        
        cluster_im(cluster_y_coord(i),cluster_x_coord(i)) = i;
        cluster_im_nor(cluster_y_coord_nor(i),cluster_x_coord_nor(i)) = i;
        
end

        
diffusion_im = diffusion_im./msd_im;
diffusion_im(isnan(diffusion_im)) = 0;

diffusion_line = diffusion_im(:);
storm_line = storm_im(:);

pixels = 1:length(storm_line);
[R,P] = corrcoef(diffusion_line,storm_line);

figure(6)
ylabels{1}='STORM Spots';
ylabels{2}='Average Diffusion Coefficient';
[ax,hlines(1),hlines(2)] = plotyy(pixels,storm_line,pixels,diffusion_line);
cfig = get(gcf,'color');
pos = [0.1  0.1  0.7  0.8];
offset = pos(3)/5.5;
hlines(1).LineWidth = 2;
hlines(2).LineWidth = 1;
hlines(1).Color = 'b';
hlines(2).Color = 'r';
set(get(ax(1),'ylabel'),'string',ylabels{1})
set(get(ax(2),'ylabel'),'string',ylabels{2})
set(ax,{'ycolor'},{'b';'r'})

title(strcat([ ' Correlation Coefficient = ', num2str(R(1,2), '%.3f'),', P-value = ', num2str(P(1,2), '%.3E')]),'FontSize',26)
xlabel('Pixel')
if repositionFigures == 1
    set(gcf,'position',[-1919,342,1920,878])
elseif repositionFigures == 2
    set(gcf,'position',[-479,1074,1901,774])
end
if saveFigures == 1
    set(gcf,'PaperPositionMode','auto')
    file1 = strcat([saveFile,'/spots_v_diff']);
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([saveFile,'/spots_v_diff.fig']);
    savefig(gcf,file1_fig)
end

storm_line_filtered = storm_line;
storm_line_filtered(storm_line_filtered>500) = 500;
[R,P] = corrcoef(storm_line_filtered,diffusion_line);

figure(7)
ylabels{1}='Filtered STORM Spots';
ylabels{2}='Average Diffusion Coefficient';
[ax,hlines(1),hlines(2)] = plotyy(pixels,storm_line_filtered,pixels,diffusion_line);
cfig = get(gcf,'color');
pos = [0.1  0.1  0.7  0.8];
offset = pos(3)/5.5;
hlines(1).LineWidth = 2;
hlines(2).LineWidth = 1;
hlines(1).Color = 'b';
hlines(2).Color = 'r';
set(get(ax(1),'ylabel'),'string',ylabels{1})
set(get(ax(2),'ylabel'),'string',ylabels{2})
set(ax,{'ycolor'},{'b';'r'})

title(strcat([ 'Filtered, Correlation Coefficient = ', num2str(R(1,2), '%.3f'),', P-value = ', num2str(P(1,2), '%.3E')]),'FontSize',26)
xlabel('Pixel')
if repositionFigures == 1
    set(gcf,'position',[-1919,523,1920,878])
elseif repositionFigures == 2
    set(gcf,'position',[-479,906,1901,774])
end
if saveFigures == 1
    set(gcf,'PaperPositionMode','auto')
    file1 = strcat([saveFile,'/spots_v_diff_filtered']);
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([saveFile,'/spots_v_diff_filtered.fig']);
    savefig(gcf,file1_fig)
end


figure(8);image(storm_im,'CDataMapping','scaled');colorbar;caxis([0 250])
title('STORM spots, Registered')
xlabel('Cell X Axis')
ylabel('Cell Y Axis')
if repositionFigures == 1
    set(gcf,'position',[-1983,331,560,420])
elseif repositionFigures == 2
    set(gcf,'position',[-479,1229,560,420])
end
if saveFigures == 1
    set(gcf,'PaperPositionMode','auto')
    file1 = strcat([saveFile,'/storm_spots_register']);
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([saveFile,'/storm_spots_register.fig']);
    savefig(gcf,file1_fig)
end

figure(9);image(msd_im,'CDataMapping','scaled');colorbar;caxis([0 6])
title('PALM spots, count')
if repositionFigures == 1
    set(gcf,'position',[-1469,331,560,420])
elseif repositionFigures == 2
    set(gcf,'position',[34,1229,560,420])
end
if saveFigures == 1
    set(gcf,'PaperPositionMode','auto')
    file1 = strcat([saveFile,'/msd_spots']);
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([saveFile,'/msd_spots.fig']);
    savefig(gcf,file1_fig)
end

figure(10);image(storm_im_nor,'CDataMapping','scaled');colorbar;caxis([0 250])
title('STORM spots, Original')
xlabel('Cell X Axis')
ylabel('Cell Y Axis')
if repositionFigures == 1
    set(gcf,'position',[-918,331,560,420])
elseif repositionFigures == 2
    set(gcf,'position',[539,1229,560,420])
end
if saveFigures == 1
    set(gcf,'PaperPositionMode','auto')
    file1 = strcat([saveFile,'/storm_spots_original']);
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([saveFile,'/storm_spots_original.fig']);
    savefig(gcf,file1_fig)
end

figure(11);image(diffusion_im./msd_im,'CDataMapping','scaled');colorbar
title('PALM spots, Diffusion Coefficients')
if repositionFigures == 1
    set(gcf,'position',[-560,331,560,420])
elseif repositionFigures == 2
    set(gcf,'position',[878,1229,560,420])
end
if saveFigures == 1
    set(gcf,'PaperPositionMode','auto')
    file1 = strcat([saveFile,'/diffusion_coefs']);
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([saveFile,'/diffusion_coefs.fig']);
    savefig(gcf,file1_fig)
end

clusterBinary = im2bw(cluster_im,1);
diffusion_in_clusters = clusterBinary.*diffusion_im;
msd_spots_in_clusters = clusterBinary .* msd_im;
not_clusters = imcomplement(clusterBinary);
diffusion_in_not_clusters = not_clusters.*diffusion_im;

clusterDiffusion = zeros(2,max(cluster_im(:))+1);
diffCoefs_in_notClusters = diffusion_in_not_clusters(:)./msd_im(:);
diffCoefs_in_notClusters(isnan(diffCoefs_in_notClusters)) = [];

clusterDiffusion(1,1) = mean(diffCoefs_in_notClusters);
clusterDiffusion(2,1) = 1;


for i = 1:max(cluster_im(:))
    if sum(cluster_im(:) == i) == 0;
        continue
    end
    
    clusterID = zeros(size(cluster_im));
    center_y = cluster_y_coord(i);
    center_x = cluster_x_coord(i);
    
    cluster_size = int32(cluster_rad(i));
    
    center_x_array = center_x-floor(cluster_size/2):center_x+floor(cluster_size/2);
    center_y_array = center_y-floor(cluster_size/2):center_y+floor(cluster_size/2);
    center_x_array(center_x_array<1) = 1;center_y_array(center_y_array<1) = 1;center_x_array(center_x_array>X) = X;center_y_array(center_y_array>Y) = Y;
    
    for j = 1:max([cluster_size,1])
        for k = 1:max([cluster_size,1])
            clusterID(center_y_array(j),center_x_array(k)) = 1;
        end
    end
            
    
    clusterDiffusion(1,i+1) = sum(sum(clusterID.*diffusion_im));
    clusterDiffusion(2,i+1) = sum(sum(clusterID.*msd_im));
    
end

clusterDiffusion(1,:) = clusterDiffusion(1,:)./clusterDiffusion(2,:);

clusterDiffusion(isnan(clusterDiffusion)) = 0;

figure(12)
barFontSize = 15;
for b = 1 : max(cluster_im(:))+1
	% Plot one single bar as a separate bar series.
	handleToThisBarSeries(b) = bar(b,clusterDiffusion(1,b)/2, 'BarWidth', 0.9);
	% Apply the color to this bar series.
    if b == 1
        set(handleToThisBarSeries(b), 'FaceColor', [1 0 0]);
    else
        set(handleToThisBarSeries(b), 'FaceColor', [0 1 0]);
    end
	hold on;
end
title('Diffusion Coefficients, in and out of Clusters')
xlabel('Cluster Number')
ylabel('Average Diffusion Coefficient')
set(gca,'XTickLabel',0:max(cluster_im(:)))

if repositionFigures == 1
    set(gcf,'position',[-560,331,560,420])
elseif repositionFigures == 2
    set(gcf,'position',[857,879,560,420])
end
if saveFigures == 1
    set(gcf,'PaperPositionMode','auto')
    file1 = strcat([saveFile,'/diff_bar']);
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([saveFile,'/diff_bar.fig']);
    savefig(gcf,file1_fig)
end

