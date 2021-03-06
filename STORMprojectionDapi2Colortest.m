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

pixelscaling = 130; % Sometimes Seongjin's data is in nm, sometimes in um. If nm, set to 130. If um, set to .130
% set "dataset" to STORM output .txt file
dataset='\\desktop-lkqeepp\G\Seongjin\AnalysisComputer2\2018_10_30_fixed_SP98_CM_Rif_KSG\SP98_NT_rRNA_Hfq\A647_5\output_fr.txt';                                          % Name of the datafile that you input....Lot of output files will have this "dataset" in their names
%set "dic_file" to image containg the dic image

dic_file = '/Users/reyer/Documents/MATLAB/SOURCE_CODES/sample_images_matt/Matt_Microscope/SP_2Colortest/DIC_5_250.tif';
dapi_file = '/Users/reyer/Documents/MATLAB/SOURCE_CODES/sample_images_matt/Matt_Microscope/SP_2Colortest/sample';
stack_dapi = imFormat_wholeImage(dapi_file,1);
dapi_objects = dapiMask(stack_dapi);

M=textread(dataset); 
red_set = 1; %Set "proper_set" to the color code you want to compute
green_set = 3;

x_col = 1;
y_col = 2;
z_col = 3;


red_particles = M(M(:,4)==red_set,:);
red_x=red_particles(:,x_col);                          % Select the column containing X coordinates in Pixel form
red_y=red_particles(:,y_col);                          % Select the column containing Y coordinates in Pixel form
red_z=red_particles(:,z_col);                          % Select the column containing Z coordinates in Pixel form

green_particles = M(M(:,4)==green_set,:);
green_x=green_particles(:,x_col);                          % Select the column containing X coordinates in Pixel form
green_y=green_particles(:,y_col);                          % Select the column containing Y coordinates in Pixel form
green_z=green_particles(:,z_col);                          % Select the column containing Z coordinates in Pixel form
 
green_start = length(red_particles)+1;

%*********Scaling the X Y Coordinates as per the pixel of the input image
filepath_dic = dic_file;
dic_image = imread(filepath_dic);                  % Put the reference image based on which you want to select the points
X = size(dic_image,1) ;                              % CapitalX is the number of X pixels in the reference image
Y = size(dic_image,2) ;                              % CapitalY is the number of Y pixels in the reference image
%scalex=X/256;                                    % Getting the scaling factors to match X and Y pixels from the image to the ones from the coordinate file
%scaley=Y/243;                                    % We are assuming that the reference image is 256 by 256 in pixel
scalex=1;
scaley=1;

red_x=(scalex.*red_x)./pixelscaling;                     % Scale it up...Pixel Scaling is to convert coordinates in nm,Angstrom to Pixel Coordinates
red_y=(scaley.*red_y)./pixelscaling;                     % Scale it up...Pixel Scaling is to convert coordinates in nm,Angstrom to Pixel Coordinates
red_z = red_z./pixelscaling; 
green_x=(scalex.*green_x)./pixelscaling;                     % Scale it up...Pixel Scaling is to convert coordinates in nm,Angstrom to Pixel Coordinates
green_y=(scaley.*green_y)./pixelscaling;                     % Scale it up...Pixel Scaling is to convert coordinates in nm,Angstrom to Pixel Coordinates
green_z = green_z./pixelscaling;


[mask,area,ellipticity,centers] = track_mask_SPtest(dic_file,pixelscaling);


% X-translation : - = Down, + = Up
% Y-translation : - = right, + = left
close all
masklg =imresize(mask,2);
imshow(masklg);hold on;scatter(red_y*2,red_x*2,0.5,'filled','r');hold on;scatter(green_y*2,green_x*2,0.5,'g');
%imshow(mask);hold on;scatter(red_y,red_x,10,'filled','r');hold on;scatter(green_y,green_x,1,'+g');
continuetranslation=1;
while continuetranslation==1
prompt={'Y_Translation(- = Down, + = Up):','X_Translation(- = right, + = left):','Y_Stretch (0-1 : Compress, >1 = Stretch','X_Stretch (0-1 : Compress, >1 = Stretch','Continue the translation(1 == Yes && 0== NO) ?'};   % A box will take in the values for the X/Ytranslation
mask_title='Translation';                             % The title of the box
answer=inputdlg(prompt,mask_title);
Xtranslation = str2num(answer{1}); 
Ytranslation = str2num(answer{2});
YStretch = str2num(answer{3});
XStretch = str2num(answer{4});
if isempty(answer{3})
    YStretch = 1;
end
if isempty(answer{4})
    XStretch = 1;
end
continuetranslation = str2num(answer{5});
red_x=(YStretch*red_x)-Xtranslation;                                % Translated
red_y=(XStretch*red_y)-Ytranslation;
green_x=(YStretch*green_x)-Xtranslation;                                % Translated
green_y=(XStretch*green_y)-Ytranslation;
imshow(masklg);hold on;scatter(red_y*2,red_x*2,0.5,'filled','r');hold on;scatter(green_y*2,green_x*2,0.5,'g');
%imshow(mask);hold on;scatter(red_y,red_x,10,'filled','r');hold on;scatter(green_y,green_x,1,'+g');
end
red_x_coord = int32(red_x);
red_y_coord = int32(red_y);
red_x_coord(red_x_coord<1) = 1;
red_y_coord(red_y_coord<1) = 1;
red_x_coord(red_x_coord>=X) = X;
red_y_coord(red_y_coord>=Y) = Y;

green_x_coord = int32(green_x);
green_y_coord = int32(green_y);
green_x_coord(green_x_coord<1) = 1;
green_y_coord(green_y_coord<1) = 1;
green_x_coord(green_x_coord>=X) = X;
green_y_coord(green_y_coord>=Y) = Y;

red_spot_coords = [(1:length(red_x_coord))' double(red_x_coord) double(red_y_coord) red_z];
red_spot_coords_fixed = red_spot_coords;

green_spot_coords = [(1:length(green_x_coord))' double(green_x_coord) double(green_y_coord) green_z];
green_spot_coords_fixed = green_spot_coords;

m2 = im2bw(mask,.01);
d2 = im2bw(dapi_objects,.01);
C = imfuse(m2,d2,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]); 
close all
C_lg =imresize(C,2);
imshow(C_lg)

continuetranslation=1;
while continuetranslation==1
prompt={'X_Translation(- = Left, + = Right):','Y_Translation(- = Down, + = Up):','X_Stretch (0-1 : Compress, >1 = Stretch','Y_Stretch (0-1 : Compress, >1 = Stretch','Continue the translation(1 == Yes && 0== NO) ?'};   % A box will take in the values for the X/Ytranslation
mask_title='Dapi (green) Translation';                             % The title of the box
answer=inputdlg(prompt,mask_title);
Xtranslation = str2num(answer{1}); 
Ytranslation = str2num(answer{2});
YStretch = round(str2num(answer{4})*Y)
XStretch = round(str2num(answer{3})*X)
if isempty(answer{4}) || str2num(answer{4}) == 0
    YStretch = size(d2,1);
end
if isempty(answer{3}) || str2num(answer{3}) == 0
    XStretch = size(d2,2);
end
continuetranslation = str2num(answer{5});
size(d2)
figure;
imshow(d2)

d2 = imtranslate(d2,[Xtranslation,-Ytranslation]);
figure;
imshow(d2)

d2 = imresize(d2,[XStretch YStretch]);
C = imfuse(m2,d2,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
C_lg=imresize(C,2);
figure;
imshow(C_lg)
end

d2 = bwlabel(C(:,:,2));
d3 = d2;
num_dapi = max(d2(:));

figure;

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
[new_mask,area,ellipticity,centers] = post_track_mask(mask,d3);

for j = 1:num_dapi
    temp_dapi_mask = d3==j;
    dapi_ids = [];
    dapi_ids = unique(temp_dapi_mask.*new_mask);
    dapi_ids(dapi_ids==0) = [];
    if isempty(dapi_ids)
        d4(d3 == j) = 0;
        continue
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
field6 = 'Dapi_Boundaries';
field7 = 'Transformed_Boundaries';
field8 = 'Transformed_Dapi_Boundaries';
field9 = 'Cell_Angle';
field10 = 'Red_Spots';
field11 = 'Green_Spots';
field12 = 'Num_Spots';
cell_struct = struct(field1,[],field2, [], field3, [], field4, [],field5, [],field6, [],field7,[],field8,[],field9,[],field10,[],field11,[],field12,[]);

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

for i = 1:max(new_mask(:))
    strcat(['Working on Cell ' num2str(i)])
    cell_struct(i).Cell = i;
    cell_struct(i).Center = [centers(i,2),centers(i,1)];
    cell_struct(i).Boundaries = bwboundaries(imdilate(new_mask==i,se));
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
    
    red_cell_spots = [];
    green_cell_spots = [];
    mask_i = imdilate(new_mask==i,se);
    dapi_mask_i = d4 == i;
    [cell_row, cell_col] = ind2sub(size(mask_i),find(mask_i));
    
    for g = 1:length(red_spot_coords)
        if sum(cell_row==red_spot_coords(g,2) & cell_col==red_spot_coords(g,3)) == 1
            red_cell_spots = [red_cell_spots red_spot_coords(g,1)];
            spot_struct(red_spot_coords(g,1)).Spot = [red_spot_coords(g,1), 1];
            spot_struct(red_spot_coords(g,1)).Coordinate = [red_x(red_spot_coords(g,1)),red_y(red_spot_coords(g,1)),red_z(red_spot_coords(g,1))];
            spot_struct(red_spot_coords(g,1)).Cell = i;
            spot_struct(red_spot_coords(g,1)).Region = 1;
        end
    end
    
    for gi = 1:length(green_spot_coords)-1
        if sum(cell_row==green_spot_coords(gi,2) & cell_col==green_spot_coords(gi,3)) == 1
            green_cell_spots = [green_cell_spots green_spot_coords(gi,1)];
            spot_struct(green_spot_coords(gi,1)+green_start).Spot = [green_spot_coords(gi,1), 3];
            spot_struct(green_spot_coords(gi,1)+green_start).Coordinate = [green_x(gi),green_y(gi),green_z(gi)];
            spot_struct(green_spot_coords(gi,1)+green_start).Region = 1;
        end
    end
    
    cell_struct(i).Num_Spots = [length(red_cell_spots),length(green_cell_spots)];
    red_cell_spot_coords = zeros(length(red_cell_spots),8);
    green_cell_spot_coords = zeros(length(green_cell_spots),8);
    
    qr = 1;
    for gs = red_cell_spots
        mask_i_temp = mask_i;
        
        for r = 1:3
            mask_i_temp = imerode(mask_i_temp,se);
            [cell_row_temp, cell_col_temp] = ind2sub(size(mask_i_temp),find(mask_i_temp));
            if sum(cell_row_temp==red_spot_coords_fixed(gs,2) & cell_col_temp==red_spot_coords_fixed(gs,3)) == 1
                spot_struct(gs).Region = spot_struct(gs).Region + 1;
            end
        end
        
        [cell_row_temp, cell_col_temp] = ind2sub(size(dapi_mask_i),find(dapi_mask_i));
        if sum(cell_row_temp==red_spot_coords_fixed(gs,2) & cell_col_temp==red_spot_coords_fixed(gs,3)) == 1
            spot_struct(gs).Region = 4;
        end
        
        
        spot_struct(gs).Transform_3D_Coordinate = spot_struct(gs).Coordinate;
        spot_struct(gs).Transform_3D_Coordinate(1) = spot_struct(gs).Transform_3D_Coordinate(1) - cell_struct(i).Center(2);
        spot_struct(gs).Transform_3D_Coordinate(2) = spot_struct(gs).Transform_3D_Coordinate(2) - cell_struct(i).Center(1);
        spot_row = spot_struct(gs).Transform_3D_Coordinate(1);
        spot_col = spot_struct(gs).Transform_3D_Coordinate(2);
        spot_struct(gs).Transform_3D_Coordinate(2) = (spot_col*sin(cell_struct(i).Cell_Angle)+spot_row*cos(cell_struct(i).Cell_Angle));
        spot_struct(gs).Transform_3D_Coordinate(1) = (spot_col*cos(cell_struct(i).Cell_Angle)-spot_row*sin(cell_struct(i).Cell_Angle));
        
        red_cell_spot_coords(qr,:) = [gs spot_struct(gs).Transform_3D_Coordinate(1) spot_struct(gs).Transform_3D_Coordinate(2) spot_struct(gs).Transform_3D_Coordinate(3) spot_struct(gs).Region spot_struct(gs).Coordinate(1) spot_struct(gs).Coordinate(2) spot_struct(gs).Coordinate(3)];
        qr = qr+1;
    end
    
    qg = 1;
    for gr = green_cell_spots
        mask_i_temp = mask_i;
        
        for r = 1:3
            mask_i_temp = imerode(mask_i_temp,se);
            [cell_row_temp, cell_col_temp] = ind2sub(size(mask_i_temp),find(mask_i_temp));
            if sum(cell_row_temp==green_spot_coords_fixed(gr,2) & cell_col_temp==green_spot_coords_fixed(gr,3)) == 1
                spot_struct(gr+green_start).Region = spot_struct(gr+green_start).Region + 1;
            end
        end
        
        [cell_row_temp, cell_col_temp] = ind2sub(size(dapi_mask_i),find(dapi_mask_i));
        if sum(cell_row_temp==green_spot_coords_fixed(gr,2) & cell_col_temp==green_spot_coords_fixed(gr,3)) == 1
            spot_struct(gr+green_start).Region = 4;
        end
        
        
        spot_struct(gr+green_start).Transform_3D_Coordinate = spot_struct(gr+green_start).Coordinate;
        spot_struct(gr+green_start).Transform_3D_Coordinate(1) = spot_struct(gr+green_start).Transform_3D_Coordinate(1) - cell_struct(i).Center(2);
        spot_struct(gr+green_start).Transform_3D_Coordinate(2) = spot_struct(gr+green_start).Transform_3D_Coordinate(2) - cell_struct(i).Center(1);
        spot_row = spot_struct(gr+green_start).Transform_3D_Coordinate(1);
        spot_col = spot_struct(gr+green_start).Transform_3D_Coordinate(2);
        spot_struct(gr+green_start).Transform_3D_Coordinate(2) = (spot_col*sin(cell_struct(i).Cell_Angle)+spot_row*cos(cell_struct(i).Cell_Angle));
        spot_struct(gr+green_start).Transform_3D_Coordinate(1) = (spot_col*cos(cell_struct(i).Cell_Angle)-spot_row*sin(cell_struct(i).Cell_Angle));
        
        green_cell_spot_coords(qg,:) = [gr+green_start spot_struct(gr+green_start).Transform_3D_Coordinate(1) spot_struct(gr+green_start).Transform_3D_Coordinate(2) spot_struct(gr+green_start).Transform_3D_Coordinate(3) spot_struct(gr+green_start).Region spot_struct(gr+green_start).Coordinate(1) spot_struct(gr+green_start).Coordinate(2) spot_struct(gr+green_start).Coordinate(3)];
        qg = qg+1;
    end
    
    red_cell_spot_coords_temp1 = red_cell_spot_coords;
    red_cell_spot_coords_temp2 = red_cell_spot_coords;
    green_cell_spot_coords_temp1 = green_cell_spot_coords;
    green_cell_spot_coords_temp2 = green_cell_spot_coords;
    
    big_axis = max(red_cell_spot_coords(:,3))-min(red_cell_spot_coords(:,3));
    small_axis = min(red_cell_spot_coords(:,2))-min(red_cell_spot_coords(:,2));
    
    if big_axis < small_axis
        red_cell_spot_coords_temp1(:,2) = red_cell_spot_coords_temp2(:,3);
        red_cell_spot_coords_temp1(:,3) = red_cell_spot_coords_temp2(:,2);
        green_cell_spot_coords_temp1(:,2) = green_cell_spot_coords_temp2(:,3);
        green_cell_spot_coords_temp1(:,3) = green_cell_spot_coords_temp2(:,2);
    end
    
    red_cell_spot_coords = red_cell_spot_coords_temp1;
    green_cell_spot_coords = green_cell_spot_coords_temp1;
    
    scatter(red_cell_spot_coords(:,2),red_cell_spot_coords(:,4),5,'filled','r');
    hold on
    scatter(green_cell_spot_coords(:,2),green_cell_spot_coords(:,4),5,'filled','g');
    grid on;title(strcat(['Cell ' num2str(i)]));xlabel('X Axis');ylabel('Z Axis')
    continuetranslation=1;
    discard_cell = 0;
    while continuetranslation==1 && all_cells_good == 0
        prompt={'X_Translation(- = Left, + = Right):','Y_Translation(+ = Up, - = Down):','Continue the translation(1 == Yes && 0== NO) ?','Discard Cell (1 == Yes) ?','Skip Rest of Cells(1 == Yes) ?'};   % A box will take in the values for the X/Ytranslation
        cell_title=strcat(['Cell ' num2str(i) ' Check']);                             % The title of the box
        answer=inputdlg(prompt,cell_title);
        Xtranslation = str2num(answer{1});
        Ztranslation = str2num(answer{2});
        continuetranslation = str2num(answer{3});
        discard_cell = str2num(answer{4});
        if discard_cell == 1
            break
        end
        all_cells_good = str2num(answer{5});
        if isempty(answer{5})
            all_cells_good = 0;
        end
        red_cell_spot_coords(:,2)=red_cell_spot_coords(:,2)+Xtranslation;                                % Translated
        red_cell_spot_coords(:,4)=red_cell_spot_coords(:,4)+Ztranslation;
        green_cell_spot_coords(:,2)=green_cell_spot_coords(:,2)+Xtranslation;                                % Translated
        green_cell_spot_coords(:,4)=green_cell_spot_coords(:,4)+Ztranslation;
        scatter(red_cell_spot_coords(:,2),red_cell_spot_coords(:,4),100,'filled','r');
        hold on
        scatter(green_cell_spot_coords(:,2),green_cell_spot_coords(:,4),100,'filled','g');
        grid on;title(strcat(['Cell ' num2str(i)]));xlabel('X Axis');ylabel('Z Axis')
    end
    
    if discard_cell == 1
        for gs = red_cell_spots
            
            spot_struct(gs).Cell = [];
            spot_struct(gs).Region = [];
            spot_struct(gs).Transform_3D_Coordinate = [];
            spot_struct(gs).Collapsed_2D_Coordinate = [];
            red_spot_coords(red_spot_coords(:,1)==gs,:) = [];
            
            
        end
        
        for gr = green_cell_spots
            
            spot_struct(gr+green_start).Cell = [];
            spot_struct(gr+green_start).Region = [];
            spot_struct(gr+green_start).Transform_3D_Coordinate = [];
            spot_struct(gr+green_start).Collapsed_2D_Coordinate = [];
            green_spot_coords(green_spot_coords(:,1)==gr,:) = [];
            
        end
        
            continue
    end
    
    cell_struct(i).Red_Spots = red_cell_spot_coords;
    cell_struct(i).Green_Spots = green_cell_spot_coords;
    
    for qt = 1:length(red_cell_spot_coords(:,1))
        
        spot_struct(red_cell_spot_coords(qt,1)).Transform_3D_Coordinate(2) = red_cell_spot_coords(qt,3);
        spot_struct(red_cell_spot_coords(qt,1)).Transform_3D_Coordinate(1) = red_cell_spot_coords(qt,2);
        spot_struct(red_cell_spot_coords(qt,1)).Transform_3D_Coordinate(3) = red_cell_spot_coords(qt,4);
        spot_struct(red_cell_spot_coords(qt,1)).Collapsed_2D_Coordinate = spot_struct(red_cell_spot_coords(qt,1)).Transform_3D_Coordinate;
        spot_struct(red_cell_spot_coords(qt,1)).Collapsed_2D_Coordinate(1) = sqrt((double(spot_struct(red_cell_spot_coords(qt,1)).Transform_3D_Coordinate(1)))^2+(double(spot_struct(red_cell_spot_coords(qt,1)).Transform_3D_Coordinate(3)))^2);
        spot_struct(red_cell_spot_coords(qt,1)).Collapsed_2D_Coordinate(3) = 0;
        spot_struct(red_cell_spot_coords(qt,1)).Distance2Center = spot_struct(red_cell_spot_coords(qt,1)).Coordinate;
        spot_struct(red_cell_spot_coords(qt,1)).Distance2Center(1) = abs(.5*spot_struct(red_cell_spot_coords(qt,1)).Transform_3D_Coordinate(1)/cell_struct(i).Cell_X_Axis);
        spot_struct(red_cell_spot_coords(qt,1)).Distance2Center(2) = abs(.5*spot_struct(red_cell_spot_coords(qt,1)).Transform_3D_Coordinate(2)/cell_struct(i).Cell_Y_Axis);
        spot_struct(red_cell_spot_coords(qt,1)).Distance2Center(3) = abs(.5*spot_struct(red_cell_spot_coords(qt,1)).Transform_3D_Coordinate(3)/cell_struct(i).Cell_X_Axis);
        spot_struct(red_cell_spot_coords(qt,1)).Distance2Membrane = Distance2Edge(cell_struct(i).Transformed_Boundaries,[spot_struct(red_cell_spot_coords(qt,1)).Collapsed_2D_Coordinate(2),spot_struct(red_cell_spot_coords(qt,1)).Collapsed_2D_Coordinate(1)]);
        red_spot_coords(red_spot_coords(:,1)==red_cell_spot_coords(qt,1),:) = [];
        
    end
    
    for gt = 1:length(green_cell_spot_coords(:,1))
        
        spot_struct(green_cell_spot_coords(gt,1)).Transform_3D_Coordinate(2) = green_cell_spot_coords(gt,3);
        spot_struct(green_cell_spot_coords(gt,1)).Transform_3D_Coordinate(1) = green_cell_spot_coords(gt,2);
        spot_struct(green_cell_spot_coords(gt,1)).Transform_3D_Coordinate(3) = green_cell_spot_coords(gt,4);
        spot_struct(green_cell_spot_coords(gt,1)).Collapsed_2D_Coordinate = spot_struct(green_cell_spot_coords(gt,1)).Transform_3D_Coordinate;
        spot_struct(green_cell_spot_coords(gt,1)).Collapsed_2D_Coordinate(1) = sqrt((double(spot_struct(green_cell_spot_coords(gt,1)).Transform_3D_Coordinate(1)))^2+(double(spot_struct(green_cell_spot_coords(gt,1)).Transform_3D_Coordinate(3)))^2);
        spot_struct(green_cell_spot_coords(gt,1)).Collapsed_2D_Coordinate(3) = 0;
        spot_struct(green_cell_spot_coords(gt,1)).Distance2Center = spot_struct(green_cell_spot_coords(gt,1)).Coordinate;
        spot_struct(green_cell_spot_coords(gt,1)).Distance2Center(1) = abs(.5*spot_struct(green_cell_spot_coords(gt,1)).Transform_3D_Coordinate(1)/cell_struct(i).Cell_X_Axis);
        spot_struct(green_cell_spot_coords(gt,1)).Distance2Center(2) = abs(.5*spot_struct(green_cell_spot_coords(gt,1)).Transform_3D_Coordinate(2)/cell_struct(i).Cell_Y_Axis);
        spot_struct(green_cell_spot_coords(gt,1)).Distance2Center(3) = abs(.5*spot_struct(green_cell_spot_coords(gt,1)).Transform_3D_Coordinate(3)/cell_struct(i).Cell_X_Axis);
        spot_struct(green_cell_spot_coords(gt,1)).Distance2Membrane = Distance2Edge(cell_struct(i).Transformed_Boundaries,[spot_struct(green_cell_spot_coords(gt,1)).Collapsed_2D_Coordinate(2),spot_struct(green_cell_spot_coords(gt,1)).Collapsed_2D_Coordinate(1)]);
        green_spot_coords(green_spot_coords(:,1)==green_cell_spot_coords(gt,1),:) = [];
        
    end
    
    
end
    
    


     
