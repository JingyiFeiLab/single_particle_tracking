clear all
clc


dim  = 3;%input('Number of D''s (2/3) : ');
channels = 1;
ref_channel = 1; % Change to most in-focus channel. Probably 2/green or 3/blue
ref_slice = 11;
slices2D = 2; % How many frames above and below reference frame (e.g. 4 = reference frame +/- 4 frames)
pix_size = .130; %Microns

int_thresh = .01; % Intensity Threshold
convolve_thresh = .05; % Threshold for Voxels to include in Convolved data
shape2D_thresh = 2;  % <--- Splitting threshold, you can change this
shape3D_thresh = 3;
conc = .3;
background_thresh = .15;
pix_neigh = floor((.08/pix_size)*12);% <- If you have no idea, try this
volume_thresh = 50;
slice_thresh = 2;
zangle_thresh = 1;
dist_thresh = 3;
low_pass_check = 0; % 1 = on. Change to 0 if you want to turn it off
gap_thresh = 5;

% Path to main file (i.e. channel) that you will use for segmentation
%filepath = strcat(['/Users/reyer/Documents/MATLAB/SOURCE_CODES/sample_images_matt/Matt_Microscope/September_3_2017_convert/manX_gfp_no_plasmid/t20/sample',num2str(cell_num)]);
filepath = '/Users/reyer/Data/STORM/05042019_Syn_data/Drug_Synthetic/DIP/sample1';


stack_o = imFormat_wholeImage(filepath,26);
        

se = [1 1 1; 1 1 1 ; 1 1 1]; % Structuring Element for basic Erosion and dilation


synth_set = zeros(1000000,3);
synth_mask = zeros(256,256,26);

pixelscaling = 130;

stack_o = stack_o./max(stack_o(:));

for i = 1:26;
    if sum(i==1:5) == 1 || sum(i==22:26) == 1 
        synth_mask(:,:,i) = imerode(bwareaopen(imdilate(imerode(bradley(stack_o(:,:,i),[pix_neigh,pix_neigh],int_thresh),se),se),50),se);
    end
    synth_mask(:,:,i) = bwareaopen(imdilate(imerode(bradley(stack_o(:,:,i),[pix_neigh,pix_neigh],int_thresh),se),se),50);
end

x_spots = 256*130*rand(1000000,1);
y_spots = 256*130*rand(1000000,1);
z_spots = 50*26*(rand(1000000,1)-.5);

synth_set = [x_spots,y_spots,z_spots];

synth_set((synth_set(:,3)< -650),:) = [];

x_coords = int32(synth_set(:,1)/pixelscaling);
y_coords = int32(synth_set(:,2)/pixelscaling);
z_coords = synth_set(:,3)/50;

x_coords(x_coords < 1) = 1;
y_coords(y_coords < 1) = 1;
x_coords(x_coords > 256) = 256;
y_coords(y_coords > 256) = 256;

synth_set_sliced = {};
new_synth_set = [];

z_level = -650;
for j = 1:26
    synth_set_sliced{j} = [];
    for i = length(synth_set):-1:1
        if (synth_set(i,3) > z_level && synth_set(i,3) <= z_level + 50)
            if synth_mask(int32(x_coords(i)),int32(y_coords(i)),j) == 1
                synth_set_sliced{j} = [synth_set_sliced{j}; x_coords(i), y_coords(i), z_coords(i)];
                new_synth_set = [new_synth_set; synth_set(i,1), synth_set(i,2), synth_set(i,3)];
            end
        end
    end
    z_level = z_level + 50;
end

dataset='/Users/reyer/Documents/MATLAB/SOURCE_CODES/sample_images_matt/Matt_Microscope/ericDapi2/output.txt';                                          % Name of the datafile that you input....Lot of output files will have this "dataset" in their names
%set "dic_file" to image containg the dic image

dic_file = '/Users/reyer/Data/STORM/05042019_Syn_data/Drug_Synthetic/DIP/DIC_1.tif';
dapi_file = '/Users/reyer/Data/STORM/05042019_Syn_data/Drug_Synthetic/DIP/DAPI1';
stack_dapi = imFormat_wholeImage(dapi_file,1);
dapi_objects = dapiMask(stack_dapi);

x_col = 1;
y_col = 2;
z_col = 3;


particles = new_synth_set;
x=particles(:,x_col);                          % Select the column containing X coordinates in Pixel form
y=particles(:,y_col);                          % Select the column containing Y coordinates in Pixel form
z=particles(:,z_col);                          % Select the column containing Z coordinates in Pixel form
                       % Select the column containing Average R for each cluster center

%*********Scaling the X Y Coordinates as per the pixel of the input image
filepath_dic = dic_file;
dic_image = imread(filepath_dic);                  % Put the reference image based on which you want to select the points
X = size(dic_image,1);                              % CapitalX is the number of X pixels in the reference image
Y = size(dic_image,2);                              % CapitalY is the number of Y pixels in the reference image
scalex=X/256;                                    % Getting the scaling factors to match X and Y pixels from the image to the ones from the coordinate file
scaley=Y/256;                                    % We are assuming that the reference image is 256 by 256 in pixel
x=(scalex.*x)./pixelscaling;                     % Scale it up...Pixel Scaling is to convert coordinates in nm,Angstrom to Pixel Coordinates
y=(scaley.*y)./pixelscaling;                     % Scale it up...Pixel Scaling is to convert coordinates in nm,Angstrom to Pixel Coordinates
z = z./pixelscaling;                  


[mask,area,ellipticity,centers] = track_mask(dic_file,pixelscaling);


% X-translation : - = Down, + = Up
% Y-translation : - = right, + = left
close all
imshow(mask);hold on;scatter(y,x,10,'filled','r')
continuetranslation=1;
while continuetranslation==1
prompt={'Y_Translation(- = Down, + = Up):','X_Translation(- = right, + = left):','Y_Stretch (0-1 : Compress, >1 = Stretch','X_Stretch (0-1 : Compress, >1 = Stretch','Continue the translation(1 == Yes && 0== NO) ?'};   % A box will take in the values for the X/Ytranslation
mask_title='Translation';                             % The title of the box
answer=inputdlg(prompt,mask_title);
Xtranslation = str2num(answer{1}); 
Ytranslation = str2num(answer{2});
YStretch = str2num(answer{3});
XStretch = str2num(answer{4});
if isempty(answer{1})
    Xtranslation = 0;
end
if isempty(answer{2})
    Ytranslation = 0;
end
if isempty(answer{3}) || YStretch == 0
    YStretch = 1;
end
if isempty(answer{4}) || XStretch == 0
    XStretch = 1;
end
continuetranslation = str2num(answer{5});
if isempty(answer{5})
    continuetranslation = 1;
end
x=(YStretch*x)-Xtranslation;                                % Translated
y=(XStretch*y)-Ytranslation;
imshow(mask);hold on;scatter(y,x,10,'filled','r')
end

x_coord = int32(x);
y_coord = int32(y);
x_coord(x_coord<1) = 1;
y_coord(y_coord<1) = 1;
x_coord(x_coord>=X) = X;
y_coord(y_coord>=Y) = Y;

spot_coords = [(1:length(x_coord))' double(x_coord) double(y_coord) z];
spot_coords_fixed = spot_coords;

m2 = im2bw(mask,.01);
d2 = im2bw(dapi_objects,.01);
C = imfuse(m2,d2,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]); 
close all
imshow(C)
continuetranslation=1;
while continuetranslation==1
prompt={'X_Translation(- = Left, + = Right):','Y_Translation(- = Down, + = Up):','X_Stretch (0-1 : Compress, >1 = Stretch','Y_Stretch (0-1 : Compress, >1 = Stretch','Continue the translation(1 == Yes && 0== NO) ?'};   % A box will take in the values for the X/Ytranslation
mask_title='Dapi (green) Translation';                             % The title of the box
answer=inputdlg(prompt,mask_title);
Xtranslation = str2num(answer{1}); 
Ytranslation = str2num(answer{2});
YStretch = round(str2num(answer{4})*Y);
XStretch = round(str2num(answer{3})*X);
if isempty(answer{1}) || str2num(answer{1}) == 0
    Xtranslation = 0;
end
if isempty(answer{2}) || str2num(answer{2}) == 0
    Ytranslation = 0;
end
if isempty(answer{4}) || str2num(answer{4}) == 0
    YStretch = size(d2,1);
end
if isempty(answer{3}) || str2num(answer{3}) == 0
    XStretch = size(d2,2);
end
continuetranslation = str2num(answer{5});
if isempty(answer{5}) 
    continuetranslation = 1;
end
d2 = imtranslate(d2,[Xtranslation,-Ytranslation]);
d2 = imresize(d2,[YStretch XStretch]);
C = imfuse(m2,d2,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
imshow(C)
end

xmask = size(mask,1);
ymask = size(mask,2);

d2 = bwlabel(C(:,:,2));
xd2 = size(d2,1);
yd2 = size(d2,2);
if xd2 > xmask
    d2(xmask+1:xd2,:) = [];
end

if yd2 > xmask
    d2(:,ymask+1:yd2,:) = [];
end
   
d3 = d2;
num_dapi = max(d2(:));

for j = 1:num_dapi
    temp_dapi_mask = d2==j;
    dapi_ids = [];
    dapi_ids = unique(temp_dapi_mask.*mask);
    dapi_ids(dapi_ids==0) = [];
    if isempty(dapi_ids)
        continue
    end
    
    if length(dapi_ids) > 1
        mask_temp_split = zeros(size(mask));
        for split_i = 1:length(dapi_ids)
            mask_temp_split = mask_temp_split + (mask==dapi_ids(split_i));
            mask(mask==dapi_ids(split_i)) = 0;
        end
        mask_temp_split = bwconvhull(mask_temp_split);
        mask = mask + dapi_ids(1)*mask_temp_split;
        dapi_ids = dapi_ids(1);
    end
   
    d3(d2 == j) = dapi_ids;
    
    
end

d3 = bwlabel(d3);
d4 = d3;
num_dapi = max(d3(:));
clear centers area ellipticity
[new_mask,area,ellipticity,centers] = post_track_mask(mask,d3,pixelscaling);

for j = 1:num_dapi
    temp_dapi_mask = d3==j;
    dapi_ids = [];
    dapi_ids = unique(temp_dapi_mask.*new_mask);
    dapi_ids(dapi_ids==0) = [];
    if isempty(dapi_ids)
        d4(d3 == j) = 0;
        continue
    end
    
    if length(dapi_ids) > 1
        mask_temp_split = zeros(size(mask));
        for split_i = 1:length(dapi_ids)
            mask_temp_split = mask_temp_split + (new_mask==dapi_ids(split_i));
            new_mask(new_mask==dapi_ids(split_i)) = 0;
        end
        mask_temp_split = bwconvhull(mask_temp_split);
        mask = mask + dapi_ids(1)*mask_temp_split;
        dapi_ids = dapi_ids(1);
    end
    
    d4(d3 == j) = dapi_ids;
    
end



close all

% colorbar


field1 = 'Cell'; % All Objects, single and multi, labeled
field2 = 'Center';
field3 = 'Cell_Y_Axis';
field4 = 'Cell_X_Axis';
field5 = 'Boundaries';
field6 = 'Expanded_Boundaries';
field7 = 'Constricted_Boundaries';
field8 = 'Dapi_Boundaries';
field9 = 'Transformed_Boundaries';
field10 = 'Transformed_Expanded_Boundaries';
field11 = 'Transformed_Constricted_Boundaries';
field12 = 'Transformed_Dapi_Boundaries';
field13 = 'Cell_Angle';
field14 = 'Spots';
field15 = 'Num_Spots';
cell_struct = struct(field1,[],field2, [], field3, [], field4, [],field5, [],field6, [],field7,[],field8,[],field9,[],field10,[],field11,[],field12,[],field13,[],field14,[],field15,[]);

field1 = 'Spot'; % All Objects, single and multi, labeled
field2 = 'Cell';
field3 = 'Coordinate';
field4 = 'Transform_3D_Coordinate'; % For Distance to Center
field5 = 'Collapsed_2D_Coordinate'; % For Distance to Membrane
field6 = 'Distance2Center'; % 
field7 = 'Distance2Membrane'; % 
field8 = 'Region';


spot_struct = struct(field1,[],field2, [], field3, [], field4, [], field5, [], field6, [], field7, [],field8,[]); 

se = [1 1 1 ; 1 1 1 ; 1 1 1];
counted = [];

all_cells_good = 0;

for i = 1:25
    strcat(['Working on Cell ' num2str(i)])
    
    if sum(new_mask(:)==i) == 0
        continue
    end
   
    
    cell_struct(i).Cell = i;
    cell_struct(i).Center = [centers(i,2),centers(i,1)];
    cell_struct(i).Boundaries = bwboundaries(new_mask==i);
    cell_struct(i).Expanded_Boundaries = bwboundaries(imdilate(new_mask==i,se));
    cell_struct(i).Constricted_Boundaries = bwboundaries(imerode(new_mask==i,se));
    cell_struct(i).Dapi_Boundaries = bwboundaries(d4==i);
    cell_struct(i).Transformed_Dapi_Boundaries = cell_struct(i).Dapi_Boundaries;
    cell_struct(i).Cell_Angle = -ellipticity(i,2)*(pi/180);
    cell_struct(i).Cell_Y_Axis = ellipticity(i,3);
    cell_struct(i).Cell_X_Axis = ellipticity(i,4);
    cell_struct(i).Transformed_Boundaries = cell_struct(i).Boundaries;
    
    cell_struct(i).Transformed_Boundaries{1,1}(:,1) = cell_struct(i).Transformed_Boundaries{1,1}(:,1) - cell_struct(i).Center(2);
    cell_struct(i).Transformed_Boundaries{1,1}(:,2) = cell_struct(i).Transformed_Boundaries{1,1}(:,2) - cell_struct(i).Center(1);
    row_border = cell_struct(i).Transformed_Boundaries{1,1}(:,1);
    col_border = cell_struct(i).Transformed_Boundaries{1,1}(:,2);
    
    cell_struct(i).Transformed_Boundaries{1,1}(:,1) = (col_border*sin(cell_struct(i).Cell_Angle)+row_border*cos(cell_struct(i).Cell_Angle));
    cell_struct(i).Transformed_Boundaries{1,1}(:,2) = (col_border*cos(cell_struct(i).Cell_Angle)-row_border*sin(cell_struct(i).Cell_Angle));
    
    cell_struct(i).Transformed_Constricted_Boundaries = cell_struct(i).Constricted_Boundaries;
    
    cell_struct(i).Transformed_Constricted_Boundaries{1,1}(:,1) = cell_struct(i).Transformed_Constricted_Boundaries{1,1}(:,1) - cell_struct(i).Center(2);
    cell_struct(i).Transformed_Constricted_Boundaries{1,1}(:,2) = cell_struct(i).Transformed_Constricted_Boundaries{1,1}(:,2) - cell_struct(i).Center(1);
    row_border = cell_struct(i).Transformed_Constricted_Boundaries{1,1}(:,1);
    col_border = cell_struct(i).Transformed_Constricted_Boundaries{1,1}(:,2);
    
    cell_struct(i).Transformed_Constricted_Boundaries{1,1}(:,1) = (col_border*sin(cell_struct(i).Cell_Angle)+row_border*cos(cell_struct(i).Cell_Angle));
    cell_struct(i).Transformed_Constricted_Boundaries{1,1}(:,2) = (col_border*cos(cell_struct(i).Cell_Angle)-row_border*sin(cell_struct(i).Cell_Angle));
    
    cell_struct(i).Transformed_Expanded_Boundaries = cell_struct(i).Expanded_Boundaries;
    
    cell_struct(i).Transformed_Expanded_Boundaries{1,1}(:,1) = cell_struct(i).Transformed_Expanded_Boundaries{1,1}(:,1) - cell_struct(i).Center(2);
    cell_struct(i).Transformed_Expanded_Boundaries{1,1}(:,2) = cell_struct(i).Transformed_Expanded_Boundaries{1,1}(:,2) - cell_struct(i).Center(1);
    row_border = cell_struct(i).Transformed_Expanded_Boundaries{1,1}(:,1);
    col_border = cell_struct(i).Transformed_Expanded_Boundaries{1,1}(:,2);
    
    cell_struct(i).Transformed_Expanded_Boundaries{1,1}(:,1) = (col_border*sin(cell_struct(i).Cell_Angle)+row_border*cos(cell_struct(i).Cell_Angle));
    cell_struct(i).Transformed_Expanded_Boundaries{1,1}(:,2) = (col_border*cos(cell_struct(i).Cell_Angle)-row_border*sin(cell_struct(i).Cell_Angle));
    
    if ~isempty(cell_struct(i).Dapi_Boundaries) 
        
        for di = 1:length(cell_struct(i).Dapi_Boundaries)
            
            cell_struct(i).Transformed_Dapi_Boundaries{di,1}(:,1) = cell_struct(i).Transformed_Dapi_Boundaries{di,1}(:,1) - cell_struct(i).Center(2);
            cell_struct(i).Transformed_Dapi_Boundaries{di,1}(:,2) = cell_struct(i).Transformed_Dapi_Boundaries{di,1}(:,2) - cell_struct(i).Center(1);
            dapi_row_border = cell_struct(i).Transformed_Dapi_Boundaries{di,1}(:,1);
            dapi_col_border = cell_struct(i).Transformed_Dapi_Boundaries{di,1}(:,2);
            cell_struct(i).Transformed_Dapi_Boundaries{di,1}(:,1) = (dapi_col_border*sin(cell_struct(i).Cell_Angle)+dapi_row_border*cos(cell_struct(i).Cell_Angle));
            cell_struct(i).Transformed_Dapi_Boundaries{di,1}(:,2) = (dapi_col_border*cos(cell_struct(i).Cell_Angle)-dapi_row_border*sin(cell_struct(i).Cell_Angle));
            
        end
    end
    
    cell_spots = [];
    mask_i = imdilate(new_mask==i,se);
    dapi_mask_i = d4 == i;
    [cell_row, cell_col] = ind2sub(size(mask_i),find(mask_i));
    
    for g = 1:length(spot_coords)
        if sum(cell_row==spot_coords(g,2) & cell_col==spot_coords(g,3)) == 1
            cell_spots = [cell_spots spot_coords(g,1)];
            spot_struct(spot_coords(g,1)).Spot = spot_coords(g,1);
            spot_struct(spot_coords(g,1)).Coordinate = [x(spot_coords(g,1)),y(spot_coords(g,1)),z(spot_coords(g,1))];
            spot_struct(spot_coords(g,1)).Cell = i;
            spot_struct(spot_coords(g,1)).Region = 1;
        end
    end
    
    cell_struct(i).Num_Spots = length(cell_spots);
    cell_spot_coords = zeros(length(cell_spots),8);
    
    qr = 1;
    for gs = cell_spots
        mask_i_temp = mask_i;
        
        for r = 1:3
            mask_i_temp = imerode(mask_i_temp,se);
            [cell_row_temp, cell_col_temp] = ind2sub(size(mask_i_temp),find(mask_i_temp));
            if sum(cell_row_temp==spot_coords_fixed(gs,2) & cell_col_temp==spot_coords_fixed(gs,3)) == 1
                spot_struct(gs).Region = spot_struct(gs).Region + 1;
            end
        end
        
        [cell_row_temp, cell_col_temp] = ind2sub(size(dapi_mask_i),find(dapi_mask_i));
        if sum(cell_row_temp==spot_coords_fixed(gs,2) & cell_col_temp==spot_coords_fixed(gs,3)) == 1
            spot_struct(gs).Region = 4;
        end
        
        
        spot_struct(gs).Transform_3D_Coordinate = spot_struct(gs).Coordinate;
        spot_struct(gs).Transform_3D_Coordinate(1) = spot_struct(gs).Transform_3D_Coordinate(1) - cell_struct(i).Center(2);
        spot_struct(gs).Transform_3D_Coordinate(2) = spot_struct(gs).Transform_3D_Coordinate(2) - cell_struct(i).Center(1);
        spot_row = spot_struct(gs).Transform_3D_Coordinate(1);
        spot_col = spot_struct(gs).Transform_3D_Coordinate(2);
        spot_struct(gs).Transform_3D_Coordinate(2) = (spot_col*sin(cell_struct(i).Cell_Angle)+spot_row*cos(cell_struct(i).Cell_Angle));
        spot_struct(gs).Transform_3D_Coordinate(1) = (spot_col*cos(cell_struct(i).Cell_Angle)-spot_row*sin(cell_struct(i).Cell_Angle));
        
        cell_spot_coords(qr,:) = [gs spot_struct(gs).Transform_3D_Coordinate(1) spot_struct(gs).Transform_3D_Coordinate(2) spot_struct(gs).Transform_3D_Coordinate(3) spot_struct(gs).Region spot_struct(gs).Coordinate(1) spot_struct(gs).Coordinate(2) spot_struct(gs).Coordinate(3)];
        qr = qr+1;
    end
    
    cell_spot_coords_temp1 = cell_spot_coords;
    cell_spot_coords_temp2 = cell_spot_coords;
    
    big_axis = max(cell_spot_coords(:,3))-min(cell_spot_coords(:,3));
    small_axis = min(cell_spot_coords(:,2))-min(cell_spot_coords(:,2));
    
    if big_axis < small_axis
        cell_spot_coords_temp1(:,2) = cell_spot_coords_temp2(:,3);
        cell_spot_coords_temp1(:,3) = cell_spot_coords_temp2(:,2);
    end
    
    cell_spot_coords = cell_spot_coords_temp1;
    
    scatter(cell_spot_coords(:,2),cell_spot_coords(:,4),100,'filled');grid on;title(strcat(['Cell ' num2str(i)]));xlabel('X Axis');ylabel('Z Axis')
    continuetranslation = 1;
    discard_cell = 0;
    while continuetranslation == 1 && all_cells_good == 0
        prompt={'X_Translation(- = Left, + = Right):','Y_Translation(+ = Up, - = Down):','Continue the translation(1 == Yes && 0== NO) ?','Discard Cell (1 == Yes) ?','Skip Rest of Cells(1 == Yes) ?'};   % A box will take in the values for the X/Ytranslation
        cell_title=strcat(['Cell ' num2str(i) ' Check']);                             % The title of the box
        answer=inputdlg(prompt,cell_title);
        Xtranslation = str2num(answer{1});
        Ztranslation = str2num(answer{2});
        continuetranslation = str2num(answer{3});
        discard_cell = str2num(answer{4});
        if isempty(answer{1}) || Xtranslation == 0
            Xtranslation = 0;
        end
        if isempty(answer{2}) || Ztranslation == 0
            Ztranslation = 0;
        end
        if discard_cell == 1
            break
        end
        all_cells_good = str2num(answer{5});
        if isempty(answer{5})
            all_cells_good = 0;
            
        end
        cell_spot_coords(:,2)=cell_spot_coords(:,2)+Xtranslation;                                % Translated
        cell_spot_coords(:,4)=cell_spot_coords(:,4)+Ztranslation;
        scatter(cell_spot_coords(:,2),cell_spot_coords(:,4),100,'filled');grid on;title(strcat(['Cell ' num2str(i)]));xlabel('X Axis');ylabel('Z Axis')
    end
    
    if discard_cell == 1
        for gs = cell_spots
            
            spot_struct(gs).Cell = [];
            spot_struct(gs).Region = [];
            spot_struct(gs).Transform_3D_Coordinate = [];
            spot_struct(gs).Collapsed_2D_Coordinate = [];
            spot_coords(spot_coords(:,1)==gs,:) = [];
            
            
        end
        
        continue
    end
    
    cell_struct(i).Spots = cell_spot_coords;
    
    for qt = 1:length(cell_spot_coords)
        
        spot_struct(cell_spot_coords(qt,1)).Transform_3D_Coordinate(2) = cell_spot_coords(qt,3);
        spot_struct(cell_spot_coords(qt,1)).Transform_3D_Coordinate(1) = cell_spot_coords(qt,2);
        spot_struct(cell_spot_coords(qt,1)).Transform_3D_Coordinate(3) = cell_spot_coords(qt,4);
        spot_struct(cell_spot_coords(qt,1)).Collapsed_2D_Coordinate = spot_struct(cell_spot_coords(qt,1)).Transform_3D_Coordinate;
        spot_struct(cell_spot_coords(qt,1)).Collapsed_2D_Coordinate(1) = sign(spot_struct(cell_spot_coords(qt,1)).Transform_3D_Coordinate(1))*sqrt((double(spot_struct(cell_spot_coords(qt,1)).Transform_3D_Coordinate(1)))^2+(double(spot_struct(cell_spot_coords(qt,1)).Transform_3D_Coordinate(3)))^2);
        spot_struct(cell_spot_coords(qt,1)).Collapsed_2D_Coordinate(3) = 0;
        spot_struct(cell_spot_coords(qt,1)).Distance2Center = spot_struct(cell_spot_coords(qt,1)).Coordinate;
        spot_struct(cell_spot_coords(qt,1)).Distance2Center(1) = abs(.5*spot_struct(cell_spot_coords(qt,1)).Transform_3D_Coordinate(1)/cell_struct(i).Cell_X_Axis);
        spot_struct(cell_spot_coords(qt,1)).Distance2Center(2) = abs(.5*spot_struct(cell_spot_coords(qt,1)).Transform_3D_Coordinate(2)/cell_struct(i).Cell_Y_Axis);
        spot_struct(cell_spot_coords(qt,1)).Distance2Center(3) = abs(.5*spot_struct(cell_spot_coords(qt,1)).Transform_3D_Coordinate(3)/cell_struct(i).Cell_X_Axis);
        spot_struct(cell_spot_coords(qt,1)).Distance2Membrane = Distance2Edge(cell_struct(i).Transformed_Boundaries,[spot_struct(cell_spot_coords(qt,1)).Collapsed_2D_Coordinate(2),spot_struct(cell_spot_coords(qt,1)).Collapsed_2D_Coordinate(1)]);
        spot_coords(spot_coords(:,1)==cell_spot_coords(qt,1),:) = [];
        
    end
    
    
end









