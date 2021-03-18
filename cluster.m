% This is a script which does the analysis of the data generated from DBSCAN script run in the master_cluster script 
% Please send me a message at dgvjay@illinois.edu if you have any doubts about these scripts- Digvijay
                                                        
rmdir('ClusterFiles','s')                                               % Please make an empty folder called ClusterFiles which basically stores .mat file for each cluster but there is rmdir term here to remove any previous set from a previous data set and make a new one for new dataset
mkdir('ClusterFiles');
Radius_Gyration=[];                                                     % Stores the 3D(x,y,z considered)    Radius/Size of each cluster formed
Radius_Gyration2D=[];                                                   % Stores the 2D(only x,y considered) Radius/Size of each cluster formed
cluster_centers=[];                                                     % Stores the x,y,z center of each cluster formed
clusterpointsnumber=0;

for i=1:max(class)
     Clusterx=[];
     Clustery=[];
     Clusterz=[];
     CurrentClusterNumber=i;                                            % Keeps track of current cluster number as we go through all the points and find out about their cluster number generated from DBSCAN data
     Cluster=[];
     for j=1:numel(class)
       if class(j)==CurrentClusterNumber    
       Clusterx=[Clusterx data(j,1)];
       Clustery=[Clustery data(j,2)];
       Clusterz=[Clusterz data(j,3)];
       end
     end
     vector=CurrentClusterNumber*ones(numel(Clusterx),1);               % TotalClusters is a data file which has 4 columns..first 3 is x,y,z cooordinate and the 4th one is the column containing the cluster number of each data point
     Clusterx=Clusterx';
     Clustery=Clustery';
     Clusterz=Clusterz';
     clusterpointsnumber=clusterpointsnumber+numel(Clusterx);
     TotalClusters=[TotalClusters; [Clusterx Clustery Clusterz vector]];% TotalClusters file as described above
     Cluster      =[Cluster      ; [Clusterx Clustery Clusterz]];
     meanclusterx=mean(Clusterx);
     meanclustery=mean(Clustery);
     meanclusterz=mean(Clusterz);
     Rg=0;
     Rg2D=0;
     for k=1:numel(Clusterx)
         Rg=Rg+(Clusterx(k)-meanclusterx).^2+(Clustery(k)-meanclustery).^2+(Clusterz(k)-meanclusterz).^2;   % Calculating the Radius of each cluster(both 2D and 3D form as described above
         Rg2D=Rg2D+(Clusterx(k)-meanclusterx).^2+(Clustery(k)-meanclustery).^2;
     end
     if numel(meanclusterx)>=1
     Radius_Gyration=[Radius_Gyration;[sqrt(Rg/numel(Clusterx)) numel(Clusterx)]];                          % Stores the Radius of Gyration and the number of points in the particule cluster...this is later stored as a SHEETFILE
     Radius_Gyration2D=[Radius_Gyration2D;[sqrt(Rg2D/numel(Clusterx)) numel(Clusterx)]];                    % Stores the Radius of Gyration_2D and the number of points in the particule cluster...this is later stored as a SHEETFILE
     cluster_centers=[cluster_centers;[meanclusterx meanclustery meanclusterz sqrt(Rg/numel(Clusterx)) numel(Clusterx)]];
     clusterfile = sprintf('./ClusterFiles/Cluster_%d.mat',CurrentClusterNumber);
     save(clusterfile,'Cluster');
     end
end
   
% copy_table = array2table(copy_array);copy_table.Properties.VariableNames = {'Time' 'Mean_GFP' 'GFP_Sigma' 'Mean_ptsG' 'ptsG_Sigma' 'Mean_SgrS' 'SgrS_Sigma'};
%     int_table = array2table(int_array);int_table.Properties.VariableNames = {'Time' 'Mean_Blue' 'Blue_Sigma' 'Blue_SE' 'Mean_Green' 'Green_Sigma' 'Green_SE' 'Mean_Red' 'Red_Sigma' 'Red_SE' 'Mean_Violet' 'Violet_Sigma'};
%     
%     copy_table_file = strcat([dateDir,'VThresh_BackSubCopy_numbers.csv']);
%     int_table_file = strcat([dateDir,'VThresh_BackSubIntensity.csv']);
% writetable(copy_table,copy_table_file);
%     writetable(int_table,int_table_file);
    


if Radius_Gyration~=0      % i.e if there are any clusters formed....there could be a case where all the points are noise.
sheetfile2=sprintf('%s_CSVFiles_ColorCode%d/Sheet_%d_%d_%d.csv',dataset,ColorCode,Npts,Eps,ColorCode);
sheetfile2D=sprintf('%s_CSVFiles_ColorCode%d/Sheet_2D_%d_%d_%d.csv',dataset,ColorCode,Npts,Eps,ColorCode);
datfile=sprintf('%s_Coordinates_ColorCode%d/Coordinates_%d_%d_%d.dat',dataset,ColorCode,Npts,Eps,ColorCode);
datfile_centers=sprintf('%s_Cluster_Centers_ColorCode%d_%d_%d.dat',dataset,ColorCode,Npts,Eps);
writetable(array2table(Radius_Gyration),sheetfile2) ;    %Stored the Radius of Gyration and the number of points in the particule cluster
writetable(array2table(Radius_Gyration2D),sheetfile2D) ; %Stored the Radius of Gyration2D and the number of points in the particule cluster
dlmwrite(datfile_centers,cluster_centers,' ');
dlmwrite( datfile,TotalClusters,' ');
averagesheet=[averagesheet;[Npts Eps mean(Radius_Gyration(:,[1])) mean(Radius_Gyration2D(:,[1])) mean(Radius_Gyration(:,[2])) max(class) (3*clusterpointsnumber/numel(data))]];
%averegasheetfile stores (six columns) : 
%Npts 
%Eps 
%Average_Size_of_Cluster_at_the_given_parameters 
%Average_Size_of_Cluster_(2D)_at_the_given_parameters 
%Average_Number_of_Points_in_each_Cluster_at_the_given_parameters 
%Number_of_Cluster_Formed_atthe_given_parameters 
%Fraction_of_points_which_are_in_clusters_atthe_givenparameters


end

%%PLOTTING THE POINTS ACCORDING TO THE CLUSTERS FORMED

cc=hsv(max(class));        % Give different colours to each cluster formed
for i=1:max(class)
   clusterfile=sprintf('./ClusterFiles/Cluster_%d.mat',i);
    load(clusterfile);
    for j=1:(numel(Cluster)/3);
    plot3(Cluster(j,1),Cluster(j,2),Cluster(j,3),'.','color',cc(i,:));
    if(i==1 && j==1)         
       hold
    end
    end
end
axis equal
% PLOTTING THE NOISE POINTS
%Storing the x,y,z coodrinates of the points considered as noise (results
%from the DBSCAN)
Noisex=[];
Noisey=[];
Noisez=[];
NoiseFound=0;
 for i=1:numel(class)
   if class(i)==-1
      NoiseFound=NoiseFound+1;
      Noisex=[Noisex data(i,1)];
      Noisey=[Noisey data(i,2)];
      Noisez=[Noisez data(i,3)];
    plot3(data(i,1),data(i,2),data(i,3),'k.');
     end
 end
axis equal
% ***** Calculating different types of distances between the cluster
% ***** min_distance_list is the minimum distance between the points of any two cluster
% ***** min_exclusion_distance_list is the minimum distance between the points of any two clusters but with a slightly exclusion meaning that if 2 is closest
% to 1 then 2 would be the listed as closest distance neighbour of 1 but if we go to point 2 and its closest neighbour also turns out to 1 then we
% list its 2nd closest neighbour distance in the list...simple meaning that a pair is never repeated in this list.
% ***** cluster_centers_distance_list is the distance between the cluster centers

% ***min_distance_list
min_distance_list=[];
for i=1:max(class)
    clusterfile1 = sprintf('./ClusterFiles/Cluster_%d.mat',i);
    load (clusterfile1);
    clust1=Cluster;
    dist.(sprintf('f_%d',i))=[];
    for j=1:max(class)
        if (j~=i)
         clusterfile2 = sprintf('./ClusterFiles/Cluster_%d.mat',j);
         load (clusterfile2);
         clust2=Cluster;
         mindistance=10000000;
         for a=1:(numel(clust1)/3)
              for b=1:(numel(clust2)/3)
               tempdist=sqrt((clust1(a,1)-clust2(b,1)).^2+(clust1(a,2)-clust2(b,2)).^2+(clust1(a,3)-clust2(b,3)).^2);    
               if (tempdist < mindistance)
                 mindistance=tempdist;
               end
              end
         end
         dist.(sprintf('f_%d',i))=[dist.(sprintf('f_%d',i));[i j mindistance]];
         min_distance_list=[min_distance_list;[i j mindistance]];
        end     
    end 
end

datfile_min_distance_list=sprintf('%s_Cluster_min_distance_list_betweenClusters_ColorCode%d_%d_%d.dat',dataset,ColorCode,Npts,Eps);
dlmwrite(datfile_min_distance_list,min_distance_list,' ');

% ***Min_Exclusion_Distance
min_exclusion_distance_list=[];
for i=1:max(class)
sortedlist=sortrows(dist.(sprintf('f_%d',i)),3);
check=0;
for a=1:numel(min_exclusion_distance_list)/3
    if(min_exclusion_distance_list(a,1)==min_exclusion_distance_list(1,2)&& min_exclusion_distance_list(a,2)==i)
        check=1;
    end
end
if(check==0)
min_exclusion_distance_list=[min_exclusion_distance_list;[i sortedlist(1,2) sortedlist(1,3)]];
else
min_exclusion_distance_list=[min_exclusion_distance_list;[i sortedlist(2,2) sortedlist(2,3)]];
end
end

datfile_min_exclusion_distance_list=sprintf('%s_Cluster_min_exclusion_distance_list_betweenClusters_ColorCode%d_%d_%d.dat',dataset,ColorCode,Npts,Eps);
dlmwrite(datfile_min_exclusion_distance_list,min_exclusion_distance_list,' ');

% ***cluster_Center_distance_list
cluster_centers_distance_list=[];
for i=1:numel(cluster_centers)/5
  for j=i+1:numel(cluster_centers)/5
   cluster_center_distance=sqrt((cluster_centers(j,1)-cluster_centers(i,1)).^2+(cluster_centers(j,2)-cluster_centers(i,2)).^2+(cluster_centers(j,3)-cluster_centers(i,3)).^2);
   cluster_centers_distance_list=[cluster_centers_distance_list ;[i j cluster_center_distance]];
  end
end

datfile_centers_distance_list=sprintf('%s_Cluster_centers_distance_list_betweenClusters_ColorCode%d_%d_%d.dat',dataset,ColorCode,Npts,Eps);
dlmwrite(datfile_centers_distance_list,cluster_centers_distance_list,' ');
