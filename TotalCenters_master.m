function Total_centers = TotalCenters_master(inputfile,Start_Npts,End_Npts,Increment_Npts,Start_Eps,End_Eps,Increment_Eps)
%inputfile is a .txt file which has colorcode in one columns and fullfilenames for clustercenters for the respective colorcode in the other column.
% 1  43C+aMG_med405_sgrS_A647_ptsG_A568_1_noTrail.txt_Cluster_Centers_ColorCode1_2_20.dat
% 2  43C+aMG_med405_sgrS_A647_ptsG_A568_1_noTrail.txt_Cluster_Centers_ColorCode2_2_20.dat
% 3  43C+aMG_med405_sgrS_A647_ptsG_A568_1_noTrail.txt_Cluster_Centers_ColorCode3_2_20.dat
% 4  43C+aMG_med405_sgrS_A647_ptsG_A568_1_noTrail.txt_Cluster_Centers_ColorCode4_2_20.dat


% Example of usage :
% TotalCenters('input.txt');


% OUTPUT: OUTPUT will be a file called "total_centers.dat"....please rename it,if you need to, accordingly.

Npts = Start_Npts:Increment_Npts:End_Npts;
Eps = Start_Eps:Increment_Eps:End_Eps;
filenames = {};


k = 0;

for i = Npts
    for j = Eps
        k = k+1;
        filenames{k} = strcat([inputfile]);
    end
end

Total_centers=[];
for i=1:numel(filenames);
    if exist(filenames{i})==0
        continue
    end
    data=textread(filenames{i});
    ColorCode_value=i;
    data=[data ColorCode_value*ones(numel(data)/5,1)];
    Total_centers=[Total_centers; data];

end

end
