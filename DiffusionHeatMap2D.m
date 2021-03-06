clear cell_center_x cell_center_z cell_color cell_region
close all

%% End length of each bin in microns. User-Determined
% You can also determine how many bins there should be, e.g.
% micron_bins = [0,3,5,7,10] -> gives 4 bins, binning cells from 0 to 3
% microns in length, 3-5 microns in length, 5-7, and 7-10

parentFile = '/Users/reyer/Data/STORM';

%dataSet = 'SP98_NT';
dataSet = 'SP151_aTC';
dataFile = strcat([parentFile,'/',dataSet]);
dates = {'12_06','12_13'}; 
micron_bins = [0,3,5,7,10];
samples_per_date = [7,8];
spot_cutoff = .85;

end_lengths = micron_bins/pixelscaling; % Convert lengths to pixels

for bin = 1:length(micron_bins)-1
    cells_in_bin = [];
    widths_in_bin = [];
    lengths_in_bin = [];
    
    
    cell_center_x = [];
    cell_center_y = [];
    cell_diffusion = [];
    cell_color = [];
    cell_region = [];
    
    for id = 1:length(dates)
        for ip = 0:samples_per_date(id)-1
            %spaceFile = strcat([dataFile,'/workspace_2018_',dates{id},'_SP126_CM_min3_Tr',num2str(ip),'.mat']);
            spaceFile = strcat([dataFile,'/workspace_2018_',dates{id},'_SP151_aTC_min3_Tr',num2str(ip),'.mat']);
            load(spaceFile)
            for n = 1:length(cell_struct)
                if cell_struct(n).Cell_Y_Axis <= end_lengths(bin+1) && cell_struct(n).Cell_Y_Axis > end_lengths(bin) && ~isempty(cell_struct(n).Spots)
                    widths_in_bin = [widths_in_bin cell_struct(n).Cell_X_Axis];
                end
            end
        end
    end
    
    for id = 1:length(dates)
        for ip = 0:samples_per_date(id)-1
            %spaceFile = strcat([dataFile,'/workspace_2018_',dates{id},'_SP126_CM_min3_Tr',num2str(ip),'.mat']);
            spaceFile = strcat([dataFile,'/workspace_2018_',dates{id},'_SP151_aTC_min3_Tr',num2str(ip),'.mat']);
            load(spaceFile)
            for n = 1:length(cell_struct)
                if cell_struct(n).Cell_Y_Axis <= end_lengths(bin+1) && cell_struct(n).Cell_Y_Axis > end_lengths(bin) && ~isempty(cell_struct(n).Spots)
                    cells_in_bin = [cells_in_bin n];
                    lengths_in_bin = [lengths_in_bin cell_struct(n).Cell_Y_Axis];
                    
                    if isempty(cell_struct(n).Spots)
                        continue
                    end
                    
                    length_normalization = end_lengths(bin+1)/cell_struct(n).Cell_Y_Axis;
                    width_normalization = mean(widths_in_bin)/cell_struct(n).Cell_X_Axis;
                    
                    for si = 1:length(cell_struct(n).Spots(:,1))
                        
                        cell_center_x = [cell_center_x cell_struct(n).Spots(si,2)*width_normalization];
                        cell_center_y = [cell_center_y cell_struct(n).Spots(si,3)*length_normalization];
                        cell_region = [cell_region cell_struct(n).Spots(si,5)];
                        cell_diffusion = [cell_diffusion cell_struct(n).Spots(si,6)];
                        
                        
                        [cell_region,region_sort] = sort(cell_region);
                        cell_center_x = cell_center_x(region_sort);
                        cell_center_y = cell_center_y(region_sort);
                        
                        
                    end
                end
                
            end
        end
    end
    
    if isempty(cells_in_bin)
       
        continue
    end
    
    x_range = int32(min(cell_center_x):max(cell_center_x));
    y_range = int32(min(cell_center_y):max(cell_center_y));
    
    binned_diffusion = zeros(length(y_range),length(x_range));
    binned_spots = zeros(length(y_range),length(x_range));
    for ix = 1:length(x_range)
        for iy = 1:length(y_range)
            temp_diffusion = [];
            for isp = 1:length(cell_center_x)
                if ix == length(x_range) && iy < length(y_range)
                    if cell_center_x(isp) > x_range(ix) && cell_center_y(isp) > y_range(iy) && cell_center_y(isp) <= y_range(iy+1)
                        temp_diffusion = [temp_diffusion cell_diffusion(isp)];
                    end
                elseif ix < length(x_range) && iy == length(y_range)
                    if cell_center_x(isp) > x_range(ix) && cell_center_x(isp) <= x_range(ix+1) && cell_center_y(isp) > y_range(iy) 
                        temp_diffusion = [temp_diffusion cell_diffusion(isp)];
                    end
                elseif ix == length(x_range) && iy == length(y_range)
                    if cell_center_x(isp) > x_range(ix) && cell_center_y(isp) > y_range(iy)
                        temp_diffusion = [temp_diffusion cell_diffusion(isp)];
                    end
                elseif cell_center_x(isp) > x_range(ix) && cell_center_x(isp) <= x_range(ix+1) && cell_center_y(isp) > y_range(iy) && cell_center_y(isp) <= y_range(iy+1)
                    temp_diffusion = [temp_diffusion cell_diffusion(isp)];
                end
            end
            if isempty(temp_diffusion)
                binned_diffusion(iy,ix) = 0;
            else
                binned_diffusion(iy,ix) = mean(temp_diffusion);
                binned_spots(iy,ix) = length(temp_diffusion);
            end
        end
    end
    
    figure(2*bin-1)
    smooth_binned = smooth2a(binned_diffusion,1,1);
    heatmap(smooth_binned,'','',[],'MaxColorValue',.03)
    graph_title = strcat(['Diffusion Coefficients for Cells between Lengths ', num2str(micron_bins(bin)), ' and ', num2str(micron_bins(bin+1)),' micron, n= ', num2str(length(cells_in_bin)),' Cells']);
    title(graph_title)
    xlabel('Cell X Axis')
    ylabel('Cell Y Axis')
    colorbar
    pbaspect([mean(widths_in_bin)/2 mean(lengths_in_bin)/2 mean(widths_in_bin)/2])
    set(gcf,'position',[835,883,868,667])
    set(gcf,'PaperPositionMode','auto')
    file1 = strcat([dataFile,'/cells',num2str(micron_bins(bin)), '_', num2str(micron_bins(bin+1))]);
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([dataFile,'/cells',num2str(micron_bins(bin)), '_', num2str(micron_bins(bin+1)),'.fig']);
    savefig(gcf,file1_fig)

    
    figure(2*bin)
    if median(binned_spots(:)) < 1
        heatmap(binned_spots)
    else
        heatmap(binned_spots,'','',[],'MaxColorValue',(1+spot_cutoff)*median(binned_spots(:)),'MinColorValue',(1-spot_cutoff)*median(binned_spots(:)))
    end
    graph_title = strcat(['Number of Spots for Cells between Lengths ', num2str(micron_bins(bin)), ' and ', num2str(micron_bins(bin+1)),' micron, n= ', num2str(length(cells_in_bin)),' Cells']);
    title(graph_title)
    xlabel('Cell X Axis')
    ylabel('Cell Y Axis')
    colorbar
    pbaspect([mean(widths_in_bin)/2 mean(lengths_in_bin)/2 mean(widths_in_bin)/2])
    set(gcf,'position',[835,883,868,667])
    set(gcf,'PaperPositionMode','auto')
    file1 = strcat([dataFile,'/cells',num2str(micron_bins(bin)), '_', num2str(micron_bins(bin+1)),'_spots']);
    print(file1,'-painters','-depsc','-r0')
    set(gcf,'PaperPositionMode','auto')
    print(file1,'-dpng','-r0')
    file1_fig = strcat([dataFile,'/cells',num2str(micron_bins(bin)), '_', num2str(micron_bins(bin+1)),'_spots.fig']);
    savefig(gcf,file1_fig)

%     hold on
%     [x_test,y_test] = calcEllipse(0,0,mean(widths_in_bin)/2,mean(lengths_in_bin)/2,0,360);
%     scatter(y_test,x_test,100,'+k')
    
    
end

