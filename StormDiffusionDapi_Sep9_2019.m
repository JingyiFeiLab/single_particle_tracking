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

pixelscaling = 173/2; % Sometimes Seongjin's data is in nm, sometimes in um. If nm, set to 130. If um, set to .130
% set "dataset" to STORM output .txt file
msd_file = '\\128.135.176.12\FEI_Lab\Seongjin\2020_01_24_TniQ_project\SP243_aTC\Tr3\onestep_MSD_list_uni_osd.dat';
dataset=msd_file;

%set "dic_file" to image containg the dic image
dic_file = '\\128.135.176.12\FEI_Lab\Seongjin\2020_01_24_TniQ_project\SP243_aTC\DIC_3_crop_paint2.tif';
dapi_file = '\\128.135.176.12\Fei_Lab\Seongjin\2020_01_24_TniQ_project\SP242_aTC\sample_1';
stack_dapi = imFormat_wholeImage(dapi_file,1);
%dapi_objects = dapiMask2(stack_dapi,.45);    %.6
dapi_objects = dapiMask(stack_dapi);
se = [1 1 1 ; 1 1 1 ; 1 1 1];

M=textread(dataset); 
proper_set = 1; %Set "proper_set" to the color code you want to compute


x_col = 1;
y_col = 2;
z_col = 3;
diff_col = 4;


particles = M;
x=particles(:,x_col);                          % Select the column containing X coordinates in Pixel form
y=particles(:,y_col);                          % Select the column containing Y coordinates in Pixel form
z=particles(:,z_col)*0;                          % Select the column containing Z coordinates in Pixel form
%DiffCoef=particles(:,diff_col);                  % Select the column containing Number of Spots for each cluster center
DiffCoef=sqrt(particles(:,diff_col))/0.00576;  % speed here

SpotID=particles(:,5);                        % Select the column containing Average R for each cluster center

%*********Scaling the X Y Coordinates as per the pixel of the input image
filepath_dic = dic_file;
dic_image = imread(filepath_dic);                  % Put the reference image based on which you want to select the points
X = size(dic_image,1) ;                              % CapitalX is the number of X pixels in the reference image
Y = size(dic_image,2) ;                              % CapitalY is the number of Y pixels in the reference image
%scalex=X/256;                                    % Getting the scaling factors to match X and Y pixels from the image to the ones from the coordinate file
%scaley=Y/256;                                    % We are assuming that the reference image is 256 by 256 in pixel
scalex=1;                                    % Getting the scaling factors to match X and Y pixels from the image to the ones from the coordinate file
scaley=1;

x=(scalex.*x)./pixelscaling*1000;           % Scale it up...Pixel Scaling is to convert coordinates in nm,Angstrom to Pixel Coordinates
y=(scaley.*y)./pixelscaling*1000;           % Scale it up...Pixel Scaling is to convert coordinates in nm,Angstrom to Pixel Coordinates
z = z./pixelscaling;                  


%[mask,area,ellipticity,centers] = track_mask(dic_file);
[mask,area,ellipticity,centers] = track_mask_SP_old(dic_file, pixelscaling);


% X-translation : - = Down, + = Up
% Y-translation : - = right, + = left
close all

masklg =imresize(mask,2);
%imshow(mask);hold on;scatter(y,x,10,'filled','r')
imshow(masklg);hold on;scatter(y*2,x*2,2,'filled','r')


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
    %imshow(mask);hold on;scatter(y,x,10,'filled','r')
    imshow(masklg);hold on;scatter(y*2,x*2,1.5,'filled','g')
end
x_coord = int32(x);
y_coord = int32(y);
x_coord(x_coord<1) = 1;
y_coord(y_coord<1) = 1;
x_coord(x_coord>=X) = X;
y_coord(y_coord>=Y) = Y;

spot_coords = [(1:length(x_coord))' double(x_coord) double(y_coord) z double(DiffCoef)];
spot_coords_fixed = spot_coords;

m2 = im2bw(mask,.01);
d2 = im2bw(dapi_objects,.01);
C = imfuse(m2,d2,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]); 
close all
%imshow(C)
C_lg =imresize(C,2);
imshow(C_lg)

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
    C_lg =imresize(C,2);
    imshow(C_lg)
    %imshow(C)
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
    
    if length(dapi_ids) > 1
        dapi_split = concave_split(temp_dapi_mask,1,.130,2.5);
        d2(d2 == j) = 0;
        d2 = d2 + dapi_split;
    end
        
end

d2 = bwlabel(imdilate(imerode(d2,se),se),4);
d3 = d2;

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

d3 = bwlabel(d3,4);
d4 = d3;
num_dapi = max(d3(:));
clear centers area ellipticity
[new_mask,area,ellipticity,centers] = post_track_mask(mask,d3,pixelscaling);

for j = 1:max(new_mask(:))
    temp_dapi_mask = new_mask ==j;
    dapi_ids = [];
    dapi_ids = unique(temp_dapi_mask.*d3);
    dapi_ids(dapi_ids==0) = [];
    if isempty(dapi_ids)
        %d4(d3 == j) = 0;
        continue
    end
    
    for ij = 1:length(dapi_ids)
        d4(d3 == dapi_ids(ij)) = j;
    end
    
%     if length(dapi_ids) > 1
%         mask_temp_split = zeros(size(mask));
%         for split_i = 1:length(dapi_ids)
%             mask_temp_split = mask_temp_split + (new_mask==dapi_ids(split_i));
%             new_mask(new_mask==dapi_ids(split_i)) = 0;
%         end
%         mask_temp_split = bwconvhull(mask_temp_split);
%         new_mask = new_mask + dapi_ids(1)*mask_temp_split;
%         dapi_ids = dapi_ids(1);
%     end
%     
%     d4(d3 == j) = dapi_ids;
    
end

m4 = im2bw(new_mask,.01);
d4l = im2bw(d4,.01);
C = imfuse(m4,d4l,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
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
field10 = 'Spots';
field11 = 'Num_Spots';
cell_struct = struct(field1,[],field2, [], field3, [], field4, [],field5, [],field6, [],field7,[],field8,[],field9,[],field10,[],field11,[]);

field1 = 'Spot'; % All Objects, single and multi, labeled
field2 = 'Cell';
field3 = 'Coordinate';
field4 = 'Transform_3D_Coordinate'; % For Distance to Center
field5 = 'Collapsed_2D_Coordinate'; % For Distance to Membrane
field6 = 'Distance2Center'; % 
field7 = 'Distance2Membrane'; % 
field8 = 'Region';
field9 = 'DiffusionCoefficient';


spot_struct = struct(field1,[],field2, [], field3, [], field4, [], field5, [], field6, [], field7, [],field8,[],field9,[]); 

se = [1 1 1 ; 1 1 1 ; 1 1 1];
counted = [];

all_cells_good = 0;

for i = 1:max(new_mask(:))
    strcat(['Working on Cell ' num2str(i)])
    
    if sum(new_mask(:)==i) == 0
        continue
    end
   
    cell_struct(i).Cell = i;
    cell_struct(i).Center = [centers(i,2),centers(i,1)];
    cell_struct(i).Boundaries = bwboundaries(new_mask==i);
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
    
    cell_spots = [];
    mask_i = imdilate(new_mask==i,se);
    dapi_mask_i = d4 == i;
    [cell_row, cell_col] = ind2sub(size(mask_i),find(mask_i));
    
    for g = 1:length(spot_coords)
        if sum(cell_row==spot_coords(g,2) & cell_col==spot_coords(g,3)) == 1
            cell_spots = [cell_spots spot_coords(g,1)];
            spot_struct(spot_coords(g,1)).Spot = spot_coords(g,1);
            spot_struct(spot_coords(g,1)).DiffusionCoefficient = spot_coords(g,5);
            spot_struct(spot_coords(g,1)).Coordinate = [x(spot_coords(g,1)),y(spot_coords(g,1)),z(spot_coords(g,1))];
            spot_struct(spot_coords(g,1)).Cell = i;
            spot_struct(spot_coords(g,1)).Region = 1;
        end
    end
    
    cell_struct(i).Num_Spots = length(cell_spots);
    cell_spot_coords = zeros(length(cell_spots),6);
    
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
        
        cell_spot_coords(qr,:) = [gs spot_struct(gs).Transform_3D_Coordinate(1) spot_struct(gs).Transform_3D_Coordinate(2) spot_struct(gs).Transform_3D_Coordinate(3) spot_struct(gs).Region spot_struct(gs).DiffusionCoefficient];
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
    
    scatter(cell_spot_coords(:,2),cell_spot_coords(:,3),5,'filled');grid on;title(strcat(['Cell ' num2str(i)]));xlabel('X Axis');ylabel('Y Axis')
    continuetranslation=1;
    discard_cell = 0;
    while continuetranslation==1 && all_cells_good == 0
        prompt={'X_Translation(- = Left, + = Right):','Y_Translation(+ = Up, - = Down):','Continue the translation(1 == Yes && 0== NO) ?','Discard Cell (1 == Yes) ?','Skip Rest of Cells(1 == Yes) ?'};   % A box will take in the values for the X/Ytranslation
        cell_title=strcat(['Cell ' num2str(i) ' Check']);                             % The title of the box
        answer=inputdlg(prompt,cell_title);
        Xtranslation = str2num(answer{1});
        Ytranslation = str2num(answer{2});
        continuetranslation = str2num(answer{3});
        discard_cell = str2num(answer{4});
        if discard_cell == 1
            break
        end
        all_cells_good = str2num(answer{5});
        if isempty(answer{5})
            all_cells_good = 0;
        end
        cell_spot_coords(:,2)=cell_spot_coords(:,2)+Xtranslation;                                % Translated
        cell_spot_coords(:,3)=cell_spot_coords(:,3)+Ytranslation;
        scatter(cell_spot_coords(:,2),cell_spot_coords(:,3),100,'filled');grid on;title(strcat(['Cell ' num2str(i)]));xlabel('X Axis');ylabel('Y Axis')
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
    
    for qt = 1:length(cell_spot_coords(:,1))
        
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
    
    


     
