clear all
close all

file_num = 1;
pixelscaling = 130;
ColorCode = 1;
Start_Npts = 2;
End_Npts = 2;
Increment_Npts = 1;
Start_Eps = 15;
End_Eps = 15;
Increment_Eps = 1;

output=strcat(['/Users/reyer/Documents/MATLAB/SOURCE_CODES/sample_images_matt/Matt_Microscope/Pol2_31521/C4/C4_STORM.csv']);
clustering_dataanalysis_2D(output,ColorCode,Start_Npts,End_Npts,Increment_Npts,Start_Eps,End_Eps,Increment_Eps,'off',0,0);