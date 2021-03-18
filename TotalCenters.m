function TotalCenters( inputfile)
%inputfile is a .txt file which has colorcode in one columns and fullfilenames for clustercenters for the respective colorcode in the other column.
% 1  43C+aMG_med405_sgrS_A647_ptsG_A568_1_noTrail.txt_Cluster_Centers_ColorCode1_2_20.dat
% 2  43C+aMG_med405_sgrS_A647_ptsG_A568_1_noTrail.txt_Cluster_Centers_ColorCode2_2_20.dat
% 3  43C+aMG_med405_sgrS_A647_ptsG_A568_1_noTrail.txt_Cluster_Centers_ColorCode3_2_20.dat
% 4  43C+aMG_med405_sgrS_A647_ptsG_A568_1_noTrail.txt_Cluster_Centers_ColorCode4_2_20.dat


% Example of usage :
% TotalCenters('input.txt');


% OUTPUT: OUTPUT will be a file called "total_centers.dat"....please rename it,if you need to, accordingly.

[ColorCodes filenames]=textread(inputfile,'%d     %s');

Total_centers=[];
for i=1:numel(filenames);
data=textread(filenames{i});
ColorCode_value=ColorCodes(i);
data=[data ColorCode_value*ones(numel(data)/5,1)];
Total_centers=[Total_centers; data];
end

dlmwrite('total_centers.dat',Total_centers,' ');
end

