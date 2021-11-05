function [cell_mask,area,ellipticity,center] = track_mask_SP(dic_path,pixelscaling)
%file_array = numbers of samples to analyze
%dic_path = path to dic images

pix_size = pixelscaling; %Microns

int_thresh = .0001; % Intensity Threshold: lower is geneerous (default=.0001)
shape2D_thresh = 3;  % <--- Splitting threshold, you can change this % higher is more generous   (default=5.5)  8

background_thresh = .1; % lower is generous (default=.1~.2) .1~.2
pix_neigh = 5; %floor((.08/pix_size)*12); <- If you have no idea, try this. Default=8 or 11  (5 worked well) 25, 5

low_pass_thresh = 0.01; % less is more generous (try 0.01 or 0.02 ) 0.01
%cell_area_thresh =900000.3; % less is more generous (around 1.3)
cell_area_thresh =1300000; % less is more generous (around 1.3) 2200000

% Path to main file (i.e. channel) that you will use for segmentation
%filepath_dic = strcat(['/Users/reyer/Documents/MATLAB/SOURCE_CODES/sample_images_matt/Matt_Microscope/August_24_17_convert/+SgrS/t20/dic',num2str(cell_num),'.tif']);
filepath_dic = dic_path;
stack_o = mat2gray(imread(filepath_dic));



se = [1 1 1; 1 1 1 ; 1 1 1]; % Structuring Element for basic Erosion and dilation


field2 = 'Objects'; % All Objects, single and multi, labeled
field3 = 'Cell_Labels';
field4 = 'Probability';
field5 = 'Mask'; % Single Cell Selections
field6 = 'All'; %All Objects, single and multi
field7 = 'Original'; % Original Image
field8 = 'Background';
field9 = 'Non_Single';
part1 = struct(field2, [], field3, [], field4, [], field5, [], field6, [], field7, [], field8, [], field9, []);  %Table for Part 1


stack2 = zeros(size(stack_o));
xdim = size(stack_o,1);
ydim = size(stack_o,2);

edge_cut = 10;



stack2 = anisodiff2D(stack_o,8,1/7,30,1);
I=stack_o;
[r,c] = size(I);
I2 = stack2;

%I_low_pass = low_pass(I2,.01); %.025
I_low_pass = low_pass(I2,low_pass_thresh); %.025
stack2 = stack2 - I_low_pass;
stack2 = stack2./max(max(stack2));

a = bwareaopen(imdilate(imerode(bradley(stack2,[pix_neigh,pix_neigh],int_thresh),se),se),20); % was 20 before
%a = a - bwareaopen(a,400);

b = bwlabel(a,4);
a_temp = a;

for i = 1:max(max(b))
    if  sum(sum(((b==i).*stack2)))/cellArea(b,i) < background_thresh
        a_temp(b==i) = 0;
        b(b == i) = 0;
    end
end

a = a_temp;
a(1:edge_cut,:) = 0; a(r-edge_cut+1:r,:) = 0; a(:,1:edge_cut) = 0; a(:,c-edge_cut+1:c) = 0;
BW = bwareaopen(a,10);
[~,BW] = edgeBreak(BW);
%BW = imfill(imclearborder(smallID(bwlabel(BW,4))),'holes');
BW = imfill(smallID(bwlabel(BW,4),20),'holes');



objects = bwlabel(BW,4);

num = max(objects(:));

ellipse_error = zeros(num,1);
test_ellipse = {};

for i = 1:num
    
    [ellipse1,test1] = ellipseError(objects,i);
    if isempty(ellipse1) == 1 || isempty(test1) == 1
        ellipse_error(i) = shape2D_thresh+1;
        continue
    else
        test_ellipse(i) = test1;
        ellipse_error(i) = ellipseTest(ellipse1,test1,cellArea(objects,i,pix_size),pix_size);
        
    end
end

for i = 1:num
    
    if ellipse_error(i) < shape2D_thresh
        objects(objects==i) = 0;
        object_temp = zeros(size(objects));
        for l = 1:length(test_ellipse{i})
            object_temp(test_ellipse{i}(l,1),test_ellipse{i}(l,2)) = i;
        end
        object_temp = imfill(object_temp);
        objects(object_temp == i) = i;
        objects = smallID(imfill(objects),20);
        
    end
end

uni_obs = unique(objects);
for numb = 1:length(unique(objects))-1
    objects(objects == uni_obs(numb+1)) = numb;
end
num = max(objects(:));

clear centers area ellipticity

mask = BW;
cell_labels = zeros(num,1);
q = 1; %Non-Single Cells
r = 1; % Single Cells
non_single = [];

%Single Cell Prediction
for i = 1:num
    [~,~,con_peaks] = edgeOptimize(objects,i);
    
    if ellipse_error(i) >= shape2D_thresh || con_peaks>=3
        mask(objects == i) = 0;
        cell_labels(i) = 1000*2+q;
        q = q+1;
        non_single = [non_single i];
    else
        
        cell_labels(i) = 1000 + r;
        r = r+1;
    end
end


%
% Structured Array with our Data
part1.Cell_Labels = cell_labels;
part1.Objects = objects;
part1.Probability = ellipse_error;
part1.Mask = mask;
part1.All = BW;
part1.Original = I;
%part1.Background = I_low_pass;
part1.Non_Single = non_single;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Cluster Splitting


field2 = 'Objects'; % All Objects, single and multi, labeled
field3 = 'Center';
field4 = 'Area';
field5 = 'Ellipticity';
field6 = 'Cell_Labels';
field7 = 'Probability';
field8 = 'Mask'; % Single Cell Selections
field9 = 'All'; % All Objects, single and multi

part2 = struct(field2, [], field3, [], field4, [], field5, [], field6, [], field7, [], field8, [], field9, []); %Table for Part 2

objects = part1.Objects;
objects2 = zeros(size(objects));
non_single = part1.Non_Single;


I2 = part1.All;

for i = non_single
    
    strcat(['splitting cell ', num2str(i), ' ... '])
    if i > length(part1.Probability)
        continue
    end
    
    clear edge_temp bound_temp
    split_im = concaveSplit_SP(objects,i,pix_size,shape2D_thresh);
    I2(objects == i) = 0;
    I2 = I2 + split_im;
end

objects2 = bwlabel(smallID(bwlabel(I2,4),20),4);
num = max(objects2(:));


for i = 1:num
    [~,~,con_peaks] = edgeOptimize(objects2,i);
    if con_peaks>=3
        objects2(objects2==i) = 0;
    end
end

objects2 = bwlabel(objects2,4);
num = max(objects2(:));

for i = 1:num
    if sum(sum(((objects2==i).*stack2)))/cellArea(objects2,i) < background_thresh
        BW(objects2 == i) = 0;
        objects2(objects2 == i) = 0;
        
    end
end

objects2 = bwlabel(objects2,4);
num = max(objects2(:));

clear centers area ellipticity

ellipse_error = zeros(num,1);
test_ellipse = {};

for i = 1:num
    
    [ellipse1,test1] = ellipseError(objects2,i);
    
    if isempty(ellipse1) == 1 || isempty(test1) == 1
        ellipse_error(i) = shape2D_thresh+1;
        continue
    else
        test_ellipse(i) = test1;
        ellipse_error(i) = ellipseTest(ellipse1,test1,cellArea(objects2,i,pix_size),pix_size);
        
    end
end

for i = 1:num
    
    if ellipse_error(i) < shape2D_thresh
        objects2(objects2==i) = 0;
        object_temp = zeros(size(objects2));
        for l = 1:length(test_ellipse{i})
            object_temp(test_ellipse{i}(l,1),test_ellipse{i}(l,2)) = i;
        end
        object_temp = imfill(object_temp);
        objects2(object_temp == i) = i;
        objects2 = smallID(imfill(objects2),20);
    else
        objects2(objects2==i) = 0;
    end
end

uni_obs = unique(objects2);
for numb = 1:length(unique(objects2))-1
    objects2(objects2 == uni_obs(numb+1)) = numb;
end
num = max(objects2(:));

clear centers area ellipticity

ellipse_error = zeros(num,1);
test_ellipse = {};

area=zeros(num,1);

%----------adding area threshould
for i=1:num
    area(i) = cellArea(objects2,i,pix_size);
    if area(i) < cell_area_thresh
        objects2(objects2==i) = 0;
    end        
end

uni_obs = unique(objects2);
for numb = 1:length(unique(objects2))-1
    objects2(objects2 == uni_obs(numb+1)) = numb;
end
num = max(objects2(:));

%------------------

for i = 1:num
    [ellipse1,test1] = ellipseError(objects2,i);
    if isempty(ellipse1) == 1 || isempty(test1) == 1
        ellipse_error(i) = shape2D_thresh+1;
        continue
    else
        test_ellipse(i) = test1;
        ellipse_error(i) = ellipseTest(ellipse1,test1,area(i),pix_size);
    end
end

centers=zeros(max(objects2(:)),2);

for i=1:num
    centers(i,:) = cellCenter(objects2,i);
end

ellipticity = zeros(num,4);

% Re-done Ellipticity Calculation
for i=1:num
    ellipticity(i,:) = cellEllipseSPT(objects2,i);
end

mask = smallID(I2,20);
cell_labels = zeros(num,1);
q = 1; %Non-Single Cells
r = 1; % Single Cells


%Single Cell Prediction
for i = 1:num
    if ellipse_error(i) >= shape2D_thresh %  Non - Single Cells hopefully
        ellipticity(i,5:8) = NaN;
        mask(objects2 == i) = 0;
        cell_labels(i) = 1000*2+q;
        q = q+1;
    else
        
        cell_labels(i) = 1000 + r;
        r = r+1;
    end
end


part2.Cell_Labels = cell_labels;
part2.Area = area ;
part2.Objects = objects2;
part2.Center = centers;
part2.Ellipticity = ellipticity;
part2.Probability = ellipse_error;
part2.Mask = mask;
part2.All = I2;
clear cells


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    

    
total_cells = part2.Objects;

indices = unique(total_cells(:));
index = 0;

for r = 2:length(indices)
    index = index + 1 ;
    total_cells(total_cells == indices(r)) = index;
end

cell_mask = total_cells;
area = part2.Area;
ellipticity = part2.Ellipticity;
center = part2.Center;
end