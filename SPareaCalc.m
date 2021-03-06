% Creates 2 output structures. 1) One entry for each cell containing its
% coordinates, center, boundaries, angle, and axis, and 2) One entry
% for each spot containing the cell ID it belongs to its 3D distance to
% center, and its distance to the membrane.

% This code will potentially take a while to run, depending on how many
% spots there are

% You should only have to change "pixelscaling", "dataset", "dic_file",
% "proper_set", and potentially the x and y columns for the STORM output
% (Seongjin's data is not consistent, so I have to check manually.
clear all
close all
graph_title = 'Dapi V A647 Correlation ';
parentDir = '/Users/reyer/Data/STORM/EM1238_minus/';
strain = {'EM1238'}; %MR Style
%strain = [0];
Date = {'June_26_2020_Evelyne_minus'};

samples = [12];

cell_start = 0;

field1 = 'Cell'; % All Objects, single and multi, labeled
field2 = 'Area';
field3 = 'Cell_Angle';
field4 = 'Cell_Y_Axis';
field5 = 'Cell_X_Axis';
cell_struct = struct(field1,[],field2, [], field3, [], field4, [],field5,[]);

for s = 1:length(samples)
    
    
    pixelscaling = 130; % Sometimes Seongjin's data is in nm, sometimes in um. If nm, set to 130. If um, set to .130
    % set "dataset" to STORM output .txt file
    
    sampleDir = strcat([parentDir,'sample',num2str(samples(s))]);
    if exist(sampleDir,'dir')~=7
        sfolder = strcat(['/sample',num2str(samples(s))]);
        mkdir(parentDir,sfolder)
    end
    
    imageDir=strcat(['/Users/reyer/Documents/MATLAB/SOURCE_CODES/sample_images_matt/Matt_Microscope/',Date{1},'/',strain{1},'/']);
    dic_file = strcat([imageDir,'dic',num2str(samples(s)),'.tif']);
    
    se = [1 1 1 ; 1 1 1 ; 1 1 1];
    
    
    %*********Scaling the X Y Coordinates as per the pixel of the input image
    filepath_dic = dic_file;
    dic_image = imread(filepath_dic);                  % Put the reference image based on which you want to select the points
    X = size(dic_image,1) ;                              % CapitalX is the number of X pixels in the reference image
    Y = size(dic_image,2) ;
    
    [mask, area,ellipticity,centers] = track_mask_SP(dic_file,pixelscaling);
    
    
    
    
    se = [1 1 1 ; 1 1 1 ; 1 1 1];
    counted = [];
    
    
    for i = 1:max(mask(:))
        strcat(['Working on Cell ' num2str(i)])
        cell_struct(i+cell_start).Cell = i;
        cell_struct(i+cell_start).Center = [centers(i,2),centers(i,1)];
        cell_struct(i+cell_start).Area = sum(mask(:)==i);
        cell_struct(i+cell_start).Cell_Angle = -ellipticity(i,2)*(pi/180);
        cell_struct(i+cell_start).Cell_Y_Axis = ellipticity(i,3);
        cell_struct(i+cell_start).Cell_X_Axis = ellipticity(i,4);
    end
    
    
    
    
end


saveFile = strcat([parentDir,'/cellStructure.mat']);
save(saveFile)



