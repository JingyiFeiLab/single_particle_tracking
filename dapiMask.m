function dapi_objects = dapiMask(stack_dapi)
dim  = 2;%input('Number of D''s (2/3) : ');
ref_channel = 1; % Change to most in-focus channel. Probably 2/green or 3/blue
ref_slice = 1; % Change this to the most in-focus slice
slices2D = 0; % How many frames above and below reference frame (e.g. 4 = reference frame +/- 4 frames)
pix_size = .130; %Microns

int_thresh = .000001; % Intensity Threshold
convolve_thresh = .05; % Threshold for Voxels to include in Convolved data
shape2D_thresh = 3.5;  % <--- Splitting threshold, you can change this
shape3D_thresh = 2;
conc = .3;
background_thresh = .1;
pix_neigh = 8; %floor((.08/pix_size)*12); <- If you have no idea, try this
volume_thresh = 40 * (2*slices2D + 1);
slice_thresh = slices2D + 1;
zangle_thresh = 1;
dist_thresh = 4;
low_pass_check = 1; % 1 = on. Change to 0 if you want to turn it off
gap_thresh = 5;


se = [1 1 1; 1 1 1 ; 1 1 1]; % Structuring Element for basic Erosion and dilation



edge_cut = 10;


strcat(['Working on Dapi Image'])
stack_dapi = anisodiff2D(stack_dapi,1,1/7,30,1);
I=stack_dapi;
[r,c] = size(I);
I2 = stack_dapi;
if dim == 2
    I_low_pass = low_pass(I2,.025);
else
    I_low_pass = low_pass(I2,.05);
end
stack_dapi = stack_dapi - I_low_pass;
stack_dapi = stack_dapi./max(max(stack_dapi));
a = bwareaopen(imdilate(imerode(bradley(stack_dapi,[pix_neigh,pix_neigh],int_thresh),se),se),20);
if low_pass_check == 1
    I_LP = I2 - I_low_pass;
    a = a .* imdilate(im2bw(I_LP,.3),se);
end
b = bwlabel(a,4);
a_temp = a;

for i = 1:max(max(b))
    if  sum(sum(((b==i).*stack_dapi)))/cellArea(b,i) < background_thresh
        a_temp(b==i) = 0;
        b(b == i) = 0;
    end
end

a = a_temp;
a = imclearborder(a);
a(1:edge_cut,:) = 0; a(r-edge_cut+1:r,:) = 0; a(:,1:edge_cut) = 0; a(:,c-edge_cut+1:c) = 0;
BW = bwareaopen(a,10);
[~,BW] = edgeBreak(BW);
BW = imfill(imclearborder(smallID(bwlabel(BW,4),10)),'holes');

objects = bwlabel(BW,4);
dapi_objects = objects;
%dapi_objects = imerode(objects,se);



    



