close all

folderTitle = '/Users/reyer/Data/STORM/06042014_ecoli_GFP/rne_images/';
stress = 'glucose';

bound_x = [];
bound_y = [];
expanded_bound_x = [];
expanded_bound_y = [];
constricted_bound_x = [];
constricted_bound_y = [];
spot_x = [];
spot_y = [];

for i = 1:length(cell_struct)
    if isempty(cell_struct(i).Spots)
        continue
    end

%     if cell_struct(i).Cell_X_Axis < 6
%         continue
%     end
    
    for j = 1:length(cell_struct(i).Transformed_Boundaries{1,1})
        bound_x = [bound_x; cell_struct(i).Transformed_Boundaries{1,1}(j,2)];
        bound_y = [bound_y; cell_struct(i).Transformed_Boundaries{1,1}(j,1)];
    end
    
%     for k = 1:length(cell_struct(i).Transformed_Expanded_Boundaries{1,1})
%         expanded_bound_x = [expanded_bound_x; cell_struct(i).Transformed_Expanded_Boundaries{1,1}(k,2)];
%         expanded_bound_y = [expanded_bound_y; cell_struct(i).Transformed_Expanded_Boundaries{1,1}(k,1)];
%     end
%     
%     for l = 1:length(cell_struct(i).Transformed_Constricted_Boundaries{1,1})
%         constricted_bound_x = [constricted_bound_x; cell_struct(i).Transformed_Constricted_Boundaries{1,1}(l,2)];
%         constricted_bound_y = [constricted_bound_y; cell_struct(i).Transformed_Constricted_Boundaries{1,1}(l,1)];
%     end
%     
%     for j = 1:length(cell_struct(i).Spots(:,2))
%         spot_x = [spot_x; spot_struct(cell_struct(i).Spots(j,1)).Transform_3D_Coordinate(1)];
%         spot_y = [spot_y; spot_struct(cell_struct(i).Spots(j,1)).Transform_3D_Coordinate(2)];
%     end
%     
    
end

%bound_x(bound_x < -6) = [];

[bincounts,bincenters] = hist(bound_x,30);

bincounts_neg = bincounts(bincenters<0);
bincounts_pos = bincounts(bincenters>0);
bincenters_neg = bincenters(bincenters<0);
bincenters_pos = bincenters(bincenters>0);

[pks_neg,locs_neg] = findpeaks(bincounts_neg);

pks_neg = abs(pks_neg);
[pks_neg,I_neg] = sort(pks_neg,'descend');
locs_neg = locs_neg(I_neg);

[pks_pos,locs_pos] = findpeaks(bincounts_pos);

pks_pos = abs(pks_pos);
[pks_pos,I_pos] = sort(pks_pos,'descend');
locs_pos = locs_pos(I_pos);

boundary = mean([abs(bincenters_neg(locs_neg(1))),abs(bincenters_pos(locs_pos(1)))]);



figure(1);histogram(bound_x,30);title(strcat([stress, ' Regular Cell Boundaries, Boundary = ', num2str(boundary)]),'FontSize',24)
set(gcf,'position',[835,883,868,667])
% file1 = strcat([folderTitle,'regular_cell_boundaries']);
% set(gcf,'PaperPositionMode','auto')
% print(file1,'-painters','-depsc','-r0')
% print(file1,'-painters','-dpdf','-r0')
% set(gcf,'PaperPositionMode','auto')
% print(file1,'-dpng','-r0')
% file1_fig = strcat([folderTitle,'regular_cell_boundaries.fig']);
% savefig(gcf,file1_fig)
% 
% % figure(2);histogram(expanded_bound_x);title('Expanded Cell Boundaries')
% % figure(3);histogram(constricted_bound_x);title('Constricted Cell Boundaries')
% figure(4);histogram(spot_x);title(strcat([stress, ' Pre-Cutoff Spot Distribution']),'FontSize',24)
% set(gcf,'position',[835,883,868,667])
% file1 = strcat([folderTitle,'preCut_spot_distribution']);
% set(gcf,'PaperPositionMode','auto')
% print(file1,'-painters','-depsc','-r0')
% print(file1,'-painters','-dpdf','-r0')
% set(gcf,'PaperPositionMode','auto')
% print(file1,'-dpng','-r0')
% file1_fig = strcat([folderTitle,'preCut_spot_distribution.fig']);
% savefig(gcf,file1_fig)












