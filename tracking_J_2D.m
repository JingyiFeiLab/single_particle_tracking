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

data_QD0=xlsread(fileName1);
data_QD=data_QD0(2:end,:);
data_QD(isnan(data_QD))=Inf;
%data_QD=data_QD(abs(data_QD(:,7))<threshold,:);
%data_QD(:,7)=data_QD(:,7)*factor;
data_QD=sortrows(data_QD,2);


x=data_QD(:,7)*100; 
y=data_QD(:,8)*100;

t=data_QD(:,2);
positionlist=[x y t];
param.mem=20;
param.dim=2;
param.good=15;
param.quiet=0;
maxdisp=1500;
result = track( positionlist, maxdisp, param )
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

x=result(:,1)/1000;
y=result(:,2)/1000;

frame=result(:,3);
particle=result(:,4);
d=length(particle); %total row number
ptotal=particle(d);%total particle number
figure
for k=1:ptotal
        dp=find(particle==k); %dimension of matrix of particle# K 
        dplength=length(dp);
        dpmax=max(dp);
        dpmin=min(dp);
        m(k)=k; %test k value
        x1=x(dpmin:dpmax);
        y1=y(dpmin:dpmax);
        %z1=z(dpmin:dpmax);
        frame1=frame(dpmin:dpmax);
        particle1=particle(dpmin:dpmax);
        trace_um{k}=[x1,y1,frame1,particle1];
        %eval(['sub_' num2str(k) '=[x1,y1,z1,frame1]']);% creating submatrix for particle# k sub_k 
        plot(x1,y1);
        title('3D Trace','fontsize',22);
        xlabel('X (\mum)','fontsize',22);
        %xlhand = get(gca,'xlabel');
        %set(xlhand,'string','X (um)','fontsize',22)
        ylabel('Y (\mum)','fontsize',22);
        %ylhand = get(gca,'ylabel');
        %set(ylhand,'string','Y (um)','fontsize',22)
        %zlabel('Z (\mum)','fontsize',22);
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


