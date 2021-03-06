% Creates 2 output structures. 1) One entry for each cell containing its
% coordinates, center, boundaries, angle, and axis, and 2) One entry
% for each spot containing the cell ID it belongs to its 3D distance to
% center, and its distance to the membrane.

% This code will potentially take a while to run, depending on how many
% spots there are

% You should only have to change "pixelscaling", "dataset", "dic_file",
% "proper_set", and potentially the x and y columns for the STORM output
% (Seongjin's data is not consistent, so I have to check manually.


pixelscaling = 130; % Sometimes Seongjin's data is in nm, sometimes in um. If nm, set to 130. If um, set to .130
% set "dataset" to STORM output .txt file
dataset='/Users/reyer/Documents/MATLAB/SOURCE_CODES/sample_images_matt/Matt_Microscope/SeongjinTracking/Two_color_0.txt';                                          % Name of the datafile that you input....Lot of output files will have this "dataset" in their names
%set "dic_file" to image containg the dic image
dic_file = '/Users/reyer/Documents/MATLAB/SOURCE_CODES/sample_images_matt/Matt_Microscope/SeongjinTracking/dic0.tif';
M=textread(dataset); 
proper_set = 1; %Set "proper_set" to the color code you want to compute

particles = M(M(:,4)==proper_set,:);
x=particles(:,2);                          % Select the column containing X coordinates in Pixel form
y=particles(:,1);                          % Select the column containing Y coordinates in Pixel form
z=particles(:,3);                          % Select the column containing Z coordinates in Pixel form
                       % Select the column containing Average R for each cluster center

%*********Scaling the X Y Coordinates as per the pixel of the input image
filepath_dic = dic_file;
dic_image = imread(filepath_dic);                  % Put the reference image based on which you want to select the points
X = size(dic_image,1) ;                              % CapitalX is the number of X pixels in the reference image
Y = size(dic_image,2) ;                              % CapitalY is the number of Y pixels in the reference image
scalex=X/256;                                    % Getting the scaling factors to match X and Y pixels from the image to the ones from the coordinate file
scaley=Y/256;                                    % We are assuming that the reference image is 256 by 256 in pixel
x=(scalex.*x)./pixelscaling;                     % Scale it up...Pixel Scaling is to convert coordinates in nm,Angstrom to Pixel Coordinates
y=(scaley.*y)./pixelscaling;                     % Scale it up...Pixel Scaling is to convert coordinates in nm,Angstrom to Pixel Coordinates
z = z./pixelscaling;                  


[mask,area,ellipticity,centers] = track_mask(dic_file);


% X-translation : - = Down, + = Up
% Y-translation : - = right, + = left
close all
imshow(mask);hold on;scatter(y,x,10,'filled','r')
continuetranslation=1;
while continuetranslation==1
prompt={'X_Translation(- = Down, + = Up):','Y_Translation(- = right, + = left):','Continue the translation(1 == Yes && 0== NO) ?'};   % A box will take in the values for the X/Ytranslation
title='Translation';                             % The title of the box
answer=inputdlg(prompt,title);
Xtranslation = str2num(answer{1}); 
Ytranslation = str2num(answer{2});
continuetranslation = str2num(answer{3});
x=x-Xtranslation;                                % Translated
y=y-Ytranslation;
imshow(mask);hold on;scatter(y,x,10,'filled','r')
end
x_coord = int32(x);
y_coord = int32(y);
x_coord(x_coord<1) = 1;
y_coord(y_coord<1) = 1;
x_coord(x_coord>=X) = X;
y_coord(y_coord>=Y) = Y;

close all

% colorbar


field1 = 'Cell'; % All Objects, single and multi, labeled
field2 = 'Center';
field3 = 'Cell_Y_Axis';
field4 = 'Cell_X_Axis';
field5 = 'Boundaries';
field6 = 'Transformed_Boundaries';
field7 = 'Cell_Angle';
cell_struct = struct(field1,[],field2, [], field3, [], field4, [],field5, [],field6, [],field7,[]);


for i = 1:max(mask(:))
    cell_struct(i).Cell = i;
    cell_struct(i).Center = [centers(i,2),centers(i,1)];
    cell_struct(i).Boundaries = bwboundaries(mask==i);
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
end

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

for i = 1:length(x)
    
    spot_struct(i).Spot = i;
    spot_struct(i).Coordinate = [x(i),y(i),z(i)]; %Row, Column, Z
    
    
    for j = 1:max(mask(:))
        mask_j = imdilate(mask==j,se);
        if ismember(i,counted)
            continue
        end
        spot_x = x_coord(i);
        spot_y = y_coord(i);
        if mask_j(spot_x,spot_y) > 0
           
            spot_struct(i).Region = 1;
            for r = 1:3
                mask_j = imerode(mask_j,se);
                if mask_j(spot_x,spot_y) > 0
                    spot_struct(i).Region = spot_struct(i).Region + 1;
                end
            end
                
            counted = [counted i];
            spot_struct(i).Cell = j;
            spot_struct(i).Transform_3D_Coordinate = spot_struct(i).Coordinate;
            spot_struct(i).Transform_3D_Coordinate(1) = spot_struct(i).Transform_3D_Coordinate(1) - cell_struct(j).Center(2);
            spot_struct(i).Transform_3D_Coordinate(2) = spot_struct(i).Transform_3D_Coordinate(2) - cell_struct(j).Center(1);
            row = spot_struct(i).Transform_3D_Coordinate(1);
            col = spot_struct(i).Transform_3D_Coordinate(2);
            spot_struct(i).Transform_3D_Coordinate(2) = (col*sin(cell_struct(j).Cell_Angle)+row*cos(cell_struct(j).Cell_Angle));
            spot_struct(i).Transform_3D_Coordinate(1) = (col*cos(cell_struct(j).Cell_Angle)-row*sin(cell_struct(j).Cell_Angle));
            spot_struct(i).Collapsed_2D_Coordinate = spot_struct(i).Transform_3D_Coordinate;
            spot_struct(i).Collapsed_2D_Coordinate(2) = sqrt((double(spot_struct(i).Transform_3D_Coordinate(2)))^2+(double(spot_struct(i).Transform_3D_Coordinate(3)))^2);
            spot_struct(i).Collapsed_2D_Coordinate(3) = 0;
            spot_struct(i).Distance2Center = spot_struct(i).Coordinate;
            spot_struct(i).Distance2Center(1) = abs(.5*spot_struct(i).Transform_3D_Coordinate(1)/cell_struct(j).Cell_X_Axis);
            spot_struct(i).Distance2Center(2) = abs(.5*spot_struct(i).Transform_3D_Coordinate(2)/cell_struct(j).Cell_Y_Axis);
            spot_struct(i).Distance2Center(3) = abs(.5*spot_struct(i).Transform_3D_Coordinate(3)/cell_struct(j).Cell_X_Axis);
            spot_struct(i).Distance2Membrane = Distance2Edge(cell_struct(j).Transformed_Boundaries,[spot_struct(i).Collapsed_2D_Coordinate(2),spot_struct(i).Collapsed_2D_Coordinate(1)]);
        end
    end
end
        
