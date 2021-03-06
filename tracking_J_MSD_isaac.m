% this code sucks, laozi yao tmd change it!

if ~exist('fileName1','var')|| isempty(fileName1)
    [userfilein, userdirin]=uigetfile({
        '*.xls','Data file (*.xls)';...
        '*.*','All Files (*.*)'},'Select the PALM file to process',...
        'D:\mydocument\MATLAB\PALM_tracking_AnalysisCode');
    fileName1=fullfile(userdirin,userfilein);
else
    if ~exist(fileName1,'file')
        fprintf('File not found: %s\n',fileName1);
        return;
    end
end

data_QD=load(fileName1);
%data_QD=xlsread(fileName1);
%data_QD=data_QD0(2:end, :);
data_QD(isnan(data_QD))=Inf;
%data_QD=data_QD(abs(data_QD(:,7))<threshold,:);
%data_QD(:,7)=data_QD(:,7)*factor;
data_QD=sortrows(data_QD,1);


x=data_QD(:,2); 
y=data_QD(:,3);
z=data_QD(:,4);
t=data_QD(:,1);
positionlist=[x y z t];
param.mem=10;
param.dim=3;
param.good=5;
param.quiet=0;
maxdisp=500;
result_tracking = track( positionlist, maxdisp, param )
%save('result','-ascii');
clear x y z t positionlist param maxdisp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D tracking
% x=data(:,5); 
% y=data(:,6);
% t=data(:,15);
% positionlist=[x y t];
% param.mem=20;
% param.dim=2;
% param.good=30;
% param.quiet=0;
% maxdisp=1000;
% result = track( positionlist, maxdisp, param )
% %save('result','-ascii');
% clear x y t positionlist param maxdisp

x_p=result_tracking(:,1)/1000;
y_p=result_tracking(:,2)/1000;
z_p=result_tracking(:,3)/1000;
frame=result_tracking(:,4);
particle=result_tracking(:,5);
d=length(particle); %total row number
ptotal=particle(d);%total particle number
figure
for k=1:ptotal
        dp=find(particle==k); %dimension of matrix of particle# K 
        dplength=length(dp);
        dpmax=max(dp);
        dpmin=min(dp);
        m(k)=k; %test k value
        x1=x_p(dpmin:dpmax);
        y1=y_p(dpmin:dpmax);
        z1=z_p(dpmin:dpmax);
        frame1=frame(dpmin:dpmax);
        particle1=particle(dpmin:dpmax);
        trace_um{k}=[x1,y1,z1,frame1,particle1];
        %eval(['sub_' num2str(k) '=[x1,y1,z1,frame1]']);% creating submatrix for particle# k sub_k 
        plot3(x1,y1,z1);
        title('3D Trace','fontsize',22);
        xlabel('X (\mum)','fontsize',22);
        %xlhand = get(gca,'xlabel');
        %set(xlhand,'string','X (um)','fontsize',22)
        ylabel('Y (\mum)','fontsize',22);
        %ylhand = get(gca,'ylabel');
        %set(ylhand,'string','Y (um)','fontsize',22)
        zlabel('Z (\mum)','fontsize',22);
        %zlhand = get(gca,'zlabel');
        %set(zlhand,'string','Z (um)','fontsize',22)
        hold all
        %clear dp dplength dpmax dpmin x1 y1 z1 frame1
end

hold off
set(0,'DefaultAxesFontSize',18)
clear x y z frame particle k d m x1 y1 z1 frame1 dp
%axis ij
axis equal

fid=fopen([userdirin 'result_tracking.txt'],'w');
fprintf(fid,'%f %f %f %d %d\n',result_tracking');    
fclose(fid);
%%%%%%%%%%%%%% revision
%ORIGINAL save('result_tracking','-ascii');
clear x y z t positionlist param maxdisp


%This m_file seperates particles according to their ID in the 'result_tracking', and
%calculates diffusion coeficient for each particle by fitting straight
%lines to the first 4 points of MSD/dt of the particle.
%of the tracking file and outputs the submatrices, e.g. sub_k is the
%submatrix of particle# k.
%load file
%make variable x y z frame number and particle ID
tic
x=result_tracking(:,1);
y=result_tracking(:,2);
z=result_tracking(:,3); % 0.79 is the correction factor for RI mismatch for water/oil for an 1.40 NA objective
frame=result_tracking(:,4);
particle=result_tracking(:,5);
d=length(particle); %total row number
ptotal=particle(d);%total particle number
tSteps=200;
clear Dif
clear MSD_Ind
figure
for j=1:ptotal
    dp=find(particle==j); %dimension of matrix of particle# j
    dplength=length(dp);
    dpmax=max(dp); % for particlej, the first row number
    dpmin=min(dp); % for particlej, the last row number
    m(j)=j; %test j value
    x1=x(dpmin:dpmax);
    y1=y(dpmin:dpmax);
    z1=z(dpmin:dpmax);
    frame1=frame(dpmin:dpmax);
    particle1=particle(dpmin:dpmax); % Added by Anne 1/14/2014
    TraceAll{j}=[x1,y1,z1,frame1,particle1]; % Modified by Anne 1/14/2014
    Trace_XYZ=[x1,y1,z1];
    % Anne added for range of trace
    TRI_tr = DelaunayTri(x1,y1,z1);
    [ch_tr v_tr] = convexHull(TRI_tr);
    TRI_tr2d = DelaunayTri(x1,y1);
    [ch_tr2d A_tr] = convexHull(TRI_tr2d);
   
    Trace_range(j,1)=v_tr;
    Trace_range(j,2)=A_tr;
    
	%%%%%%%%%%%% this code sucks, laozi yao tmd change it! %%%%%%%%%%%%%%%%
    [MSD00,d2r0,counts]=fMSD_vect(x1,y1,z1,frame1,dpmax,dpmin,tSteps); % Call function MSD, MSD.m file must be in the same folder as this file
	MSDall{j} = MSD00; % laozi tmd added this line of code, tmd what code is this,waste laozi so much energy to change want to die....
	
    cutoff=10;  %cutoff for data points
    ind=find(counts>cutoff-1); % find index of counts that is equate to and above cutoff
    MSDCF=MSD00(ind); % find the MSD for those index from last line
    %ind1=[0,ind]; % Add (0,0) as the first point of the curve
    %MSDCF1=[0,MSDCF]; % Add (0,0) as the first point of the curve
    % ind1=ind(1:4);
    % MSDCF1=MSDCF(1:4);
    indlength=length(ind);
    if indlength>=4
        ind1=ind(1:4);
        MSDCF1=MSDCF(1:4);
        Dif(j,:)= polyfit(ind1',MSDCF1',1);
        
% ORIGINAL        f = fit(ind1',MSDCF1','poly1');
% ORIGINAL       Dif(j,:)=coeffvalues(f); % First column is the diffusion coeff (slope). To get um^2/s multiply by 2*10^(-5), to get cm^2/s, multiply by 2*10^(-13)
    else Dif(j,:)=[0,0]; ind1=ind; MSDCF1=MSDCF; %%% Changed on 2013/05/13, Dif has to have 2 columns, and if indlenth<4, ind1 and MSDCF1 has to have values
    end
    particle2(j)=j; % Added by Anne, 1/14/2014
    plot(ind,MSDCF)
    hold all
    %eval(['MSD00_' num2str(j) '=[MSD00]']);
    %eval(['d2r0_' num2str(j) '=[d2r0]']);
    %         clear dp dplength dpmax dpmin x1 y1 z1 frame1
    inst_msd{j}= inst_MSD( x1,y1,z1,frame1);
end
hold off
saveas(gcf,strcat(userdirin,'MSD'),'fig')

dt=0.015; %Default time step 0.05s=50ms
Dif1=Dif(:,1);
Dif_positive=Dif1(Dif1>0);
Dif_time=Dif_positive/(dt*2*3); %D for time step dt, for 3D
Dif_track=Dif_positive/(dt*2*3*10^6);
Dif_track_ID=particle2(Dif1>0); % Added by Anne, 1/14/2014, for trace ID of corresponding D
Dif_track_W_ID(:,1)=Dif_track;
Dif_track_W_ID(:,2)=Dif_track_ID;
figure
hist(Dif_track);
saveas(gcf,strcat(userdirin,'Dif_track'),'fig')
save(strcat(userdirin,'Inst_diffusion_coefficients.mat'),'inst_msd');
save(strcat(userdirin,'diffusion_coefficients.mat'),'Dif_track_W_ID');
save trace.um
clear x y z frame particle j d m x1 y1 z1 frame1 dp
