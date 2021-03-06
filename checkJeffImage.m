folderTitle = '/Users/reyer/Data/STORM/';

dataset='/Users/reyer/Documents/MATLAB/SOURCE_CODES/sample_images_matt/Matt_Microscope/FISH_control_002_RNA.txt';                                          % Name of the datafile that you input....Lot of output files will have this "dataset" in their names
%set "dic_file" to image containg the dic image

M=textread(dataset); 
proper_set = 3; %Set "proper_set" to the color code you want to compute

x_col = 2;
y_col = 1;
particles = M;
x_coords=particles(:,x_col);                          % Select the column containing X coordinates in Pixel form
y_coords=particles(:,y_col);                          % Select the column containing Y coordinates in Pixel form

check = [];

bit_mult = 65536;
rand_image = zeros(1024,1024);
num_spots = length(x_coords);

for i = 1:num_spots;
    
    rand_image(1024-(x_coords(i)),(y_coords(i)+1)) = 1;
        
end

figure(1);imshow(rand_image)
