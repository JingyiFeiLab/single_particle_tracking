clear all
%close all

sample = 'C7';
saveFile = strcat(['/Users/reyer/Documents/MATLAB/SOURCE_CODES/sample_images_matt/Matt_Microscope/2021_March2nd_SC35_M3/',sample,'/']);

repositionFigures = 1;
saveFigures = 0;

gfp_path=strcat(['/Users/reyer/Data/mettl3/split_image/2021_10_13_B_002_gfp.tif']);
gfp_post_path=strcat(['/Users/reyer/Documents/MATLAB/SOURCE_CODES/sample_images_matt/Matt_Microscope/2021_March2nd_SC35_M3/',sample,'/',sample,'_M3_F6_488_POST.tif']);

low_pass_check = 1;
% Path to main file (i.e. channel) that you will use for segmentation
%filepath_dic = strcat(['/Users/reyer/Documents/MATLAB/SOURCE_CODES/sample_images_matt/Matt_Microscope/August_24_17_convert/+SgrS/t20/dic',num2str(cell_num),'.tif']);

int_thresh = .0001; % Intensity Threshold
shape2D_thresh = 2;  % <--- Splitting threshold, you can change this % higher is more generous   (default=5.5)  8

background_thresh = .2; % lower is generous (default=.1~.2) .1
pix_neigh = 5;%floor((.08/pix_size)*12);% <- If you have no idea, try this. Default=8 or 11  (5 worked well) 25

low_pass_thresh = 0.025; % less is more generous (try 0.01 or 0.02 )
high_pass_thresh = 0.075;
se = [1 1 1; 1 1 1 ; 1 1 1];

stack_o = mat2gray(imread(gfp_path));
stack_2 = mat2gray(imread(gfp_post_path));

fixed_low_pass = low_pass(stack_o,low_pass_thresh); %.025
stack_o = stack_o - fixed_low_pass;
stack_o = stack_o./max(stack_o(:));
%fixed_low_pass = fixed_low_pass./max(fixed_low_pass(:));

gfp_mask = im2bw(stack_o,.1);
figure(1);imshow(gfp_mask)

% stack_fixed = mat2gray(imread(gfp_path));
% stack_moving = mat2gray(imread(gfp_post_path));
% 
% stack2 = anisodiff2D(stack_o,1,1/7,30,1);
% I=stack_o;
% [r,c] = size(I);
% I2 = stack2;
% I_low_pass = low_pass(I2,.025);
% 
% stack2 = stack2 - I_low_pass;
% stack2 = stack2./max(max(stack2));
% a = bwareaopen(imdilate(imerode(bradley(stack2,[pix_neigh,pix_neigh],int_thresh),se),se),50);
% 
% b = bwlabel(a,4);
% a_temp = a;
% 
% for i = 1:max(max(b))
%     if  sum(sum(((b==i).*stack2)))/cellArea(b,i) < background_thresh
%         a_temp(b==i) = 0;
%         b(b == i) = 0;
%     end
% end
% 
% a = a_temp;
% BW = bwareaopen(a,10);
% 
% 
% 
% objects = bwlabel(BW,4);
% 
% X = size(stack_fixed,1) ;                              % CapitalX is the number of X pixels in the reference image
% Y = size(stack_fixed,2) ; 
% 
% se = [1 1 1; 1 1 1 ; 1 1 1]; % Structuring Element for basic Erosion and dilation
% 
% fixed_low_pass = low_pass(stack_fixed,low_pass_thresh); %.025
% fixed_low_pass = fixed_low_pass./max(fixed_low_pass(:));
% moving_low_pass = low_pass(stack_moving,low_pass_thresh); %.025
% moving_low_pass = moving_low_pass./max(moving_low_pass(:));
% 
% [optimizer,metric] = imregconfig('multimodal');
% stack_recovered = imregister(moving_low_pass,fixed_low_pass,'rigid',optimizer,metric);
% tform = imregtform(moving_low_pass,fixed_low_pass,'rigid', optimizer, metric);
% register_low_pass = imwarp(moving_low_pass,tform,'OutputView',imref2d(size(fixed_low_pass)));
% register_storm = imwarp(stack_moving,tform,'OutputView',imref2d(size(fixed_low_pass)));

