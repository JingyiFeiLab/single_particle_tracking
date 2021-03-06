function dapi_objects = dapiMask2(stack_dapi,int_thresh)
s = mat2gray(stack_dapi);
se = [1 1 1; 1 1 1 ; 1 1 1];

%dapi_objects = imerode(imdilate(imfill(im2bw(s,int_thresh),'holes'),se),se);
dapi_objects = imdilate(imerode(imfill(im2bw(s,int_thresh),'holes'),se),se);



    



