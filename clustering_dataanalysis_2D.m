function clustering_dataanalysis_2D(datafile,ColorCode,Start_Npts,End_Npts,Increment_Npts,Start_Eps,End_Eps,Increment_Eps,Window_Status,imagestatus,moviestatus)

% This is a master function which runs the DBSCAN script followed by the cluster analysis script (cluster.m) which is used to do analysis on the results from the clustering analysis on any input data
% Please send me a message at dgvjay@illinois.edu if you have any doubts about these scripts- Digvijay

%datafile    = Input Filename
%ColorCode   = 1 or 2 or 3 etc  ColorCode in the dataset that you want to do clustering on...mention 0 if you want to look at all the points in the input file
%Start_Npts  = Starting Value of Npts parameter
%End_Npts    = Ending Value of Npts parameter
%Increment_Npts= Stepwise increment from Start_Npts to End_Npts
%Start_Eps   = Starting Value of Eps parameter
%End_Eps     = Ending Value of Eps parameter
%Increment_Eps=Stepwise increment from Start_Eps to End_Eps
%Window_Status= 'on' if you want to see the figure/windows during the
%                execution of the function

%               'off' if you want to supress the figure/windows during the
%                execution of the function
%imagestatus =1 if you want to make and save images as well
%moviestatus =1 if you want to make and save movies as well

% Example of USAGE :
%clustering_dataanalysis('031612_HeLA_SC35_A568_and_SC35_A647_test_1.txt',1,5,5,1,20,20,1,'off',0,0);
%clustering_dataanalysis(datafile,ColorCode,Start_Npts,End_Npts,Increment_Npts,Start_Eps,End_Eps,Increment_Eps,Window_Status,imagestatus,moviestatus)


if exist('ClusterFiles','dir')==0
    mkdir('ClusterFiles');  
end


dataset=datafile;                                          % Name of the datafile that you input....Lot of output files will have this "dataset" in their names
M=csvread(dataset,1,1);                                       % Reads in the dataset
data=[];                                                   % So here we are assuming that x,y,z coordinates come from three different channels,colorcode1 and colorcode 3 which would be mentioned in the fourth column of input data files
for t=1:size(M,1)
    
        
    data=[data ;[M(t,1) M(t,2) 0]];                       % So here's how I have saved the points belonging to a colorcode from the original dataset
        
end

if ColorCode==0
    data=M;
end

 
rmdir('ClusterFiles','s')                                      
mkdir('ClusterFiles'); 
fold1=sprintf('%s_CSVFiles_ColorCode%d',dataset,ColorCode);    % A Folder for storing CSV Files with the dataset name
fold2=sprintf('%s_Coordinates_ColorCode%d',dataset,ColorCode); % A Folder for storing Coordinate Files with the dataset name
mkdir(fold1);
mkdir(fold2);
averagesheet=[];
finalminseparations=[];
Total_Average_Separations=[];
averagefile=sprintf('%s_ColorCode%d_AverageSheet_2D.csv',dataset,ColorCode);
for Npts=Start_Npts:Increment_Npts:End_Npts                    % The Npts parameter of the DBSCAN
    for Eps=Start_Eps:Increment_Eps:End_Eps                    % The Eps  parameter of the DBSCAN
        [class type]=dbscan(data,Npts,Eps);
        TotalClusters=[];
        imagewindow = figure('visible',Window_Status);
        run ./cluster.m                                        % A script which does all the required analysis...so any additional analysis or modifications in analysis required has to made in this script
        imagefile     = sprintf('%s_%d_%d_ColorCode%d.jpg',dataset,Npts,Eps,ColorCode);
        imagefile_fig = sprintf('%s_%d_%d_ColorCode%d.fig',dataset,Npts,Eps,ColorCode);
        moviefile     = sprintf('%s_%d_%d_ColorCode%d.avi',dataset,Npts,Eps,ColorCode);
        if imagestatus==1
        saveas(imagewindow,imagefile);
        saveas(imagewindow,imagefile_fig);
        end
        if moviestatus==1
          [az,el] = view;
          rots = [0:5:360];
          for i = 1:72
            view(rots(i),el);
            MO(i) = getframe(gcf);
          end
            movie2avi(MO, moviefile, 'compression', 'None');
        end
        close(imagewindow) ;
    end
end
if isempty(averagesheet)==0
    xlswrite(averagefile,averagesheet);
end
rmdir('ClusterFiles','s');  
end




