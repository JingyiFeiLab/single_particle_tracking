close all

parentDir = '/Users/reyer/Data/STORM/WT';
strain = {'WT'}; %MR Style
pm = '+';

samples = 2;
cell_size = [];
cell_length = [];
cell_ellipse = [];

j = 0;

for i = 1:length(All_struct)
    
%     if All_struct(i).Cell_Y_Axis<15
%         j = j+1;
%         continue
%     end
    
    cell_size = [cell_size All_struct(i).Area];
    cell_length = [cell_length All_struct(i).Cell_Y_Axis];
    cell_ellipse = [cell_ellipse ellipticity(i,1)];
    %cell_width = [cell_width All_struct(i).Cell_X_Axis];
          
end

figure(1)
histogram(cell_size)
title(strcat([strain{1},' cell size = ', num2str(mean(cell_size)), ' (', num2str(std(cell_size)), ')', ' n = ', num2str(length(All_struct)-j)]))
set(gca, 'FontName', 'Arial')
set(gca,'linewidth',1)
%set(gcf,'position',[835,883,868,667])
set(gcf,'position',[854,226,310,256])
file1 = strcat([parentDir,'/cell_size']);
set(gcf,'PaperPositionMode','auto')
saveas(gcf,file1,'png')

figure(2)
histogram(cell_length)
title(strcat([strain{1},' cell length = ', num2str(mean(cell_length)), ' (', num2str(std(cell_length)), ')', ' n = ', num2str(length(All_struct)-j)]))
set(gca, 'FontName', 'Arial')
set(gca,'linewidth',1)
%set(gcf,'position',[835,883,868,667])
set(gcf,'position',[854,226,310,256])
file1 = strcat([parentDir,'/cell_length']);
set(gcf,'PaperPositionMode','auto')
saveas(gcf,file1,'png')


figure(3)
histogram(cell_ellipse)
title(strcat([strain{1},' cell ellipticity = ', num2str(mean(cell_ellipse)), ' (', num2str(std(cell_ellipse)), ')' , ' n = ', num2str(length(All_struct)-j)]))
set(gca, 'FontName', 'Arial')
set(gca,'linewidth',1)
%set(gcf,'position',[835,883,868,667])
set(gcf,'position',[854,226,310,256])
file1 = strcat([parentDir,'/cell_width']);
saveas(gcf,file1,'png')