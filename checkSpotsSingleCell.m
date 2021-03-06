close all

parentDir = '/Users/reyer/Data/STORM/';

all_cell_ids = [];
perma_cell_ids = [];

strain = 156; %MR Style
%strain = [0];
Date = 'September_7_2019';
%Date = {'eric_images'};
time = 30;
sample = 0;
color = 'red';

cell_check = 31;


for i = 1:length(cell_region_struct)
    all_cell_ids = [all_cell_ids cell_region_struct(i).Cell];
    perma_cell_ids = [perma_cell_ids cell_region_struct(i).Cell];
end

for j = 1
    
    one_cell_id = cell_check;
    
    
    close all
    figure(1)
    
    
    %one_cell_id = 10; % If you want to create a scatter plot for just one cell, put it's ID here
    plot_all_cells = 0; % If you want to create a heatmap for all spots in all cells, set to 1
    percent_from_pole = 0; % Cutoff spots within this percentage of cell axis of cell pole
    
    cell_center_x = [];
    cell_center_y = [];
    cell_center_x2d = [];
    cell_center_y2d = [];
    cell_center_z = [];
    cell_color = [];
    cell_region = [];
    
    from_pole_temp = [];
    
    levels = 10; %How many heights for cell spot density map? More takes longer
    % 10 takes a little bit of time, 100 a lot
    
    for si = 1:length(cell_struct(one_cell_id).Spots(:,1))
        spot_id = cell_struct(one_cell_id).Spots(si,1);
        spot_row = spot_struct(spot_id).Collapsed_2D_Coordinate(1);
        spot_col = spot_struct(spot_id).Collapsed_2D_Coordinate(2);
        new_spot_row = int32((spot_row*sin(-1*cell_struct(one_cell_id).Cell_Angle)+spot_col*cos(-1*cell_struct(one_cell_id).Cell_Angle)) + cell_struct(one_cell_id).Center(2));
        new_spot_col = int32((spot_row*cos(-1*cell_struct(one_cell_id).Cell_Angle)-spot_col*sin(-1*cell_struct(one_cell_id).Cell_Angle)) + cell_struct(one_cell_id).Center(1));
        
        if cell_struct(one_cell_id).Spots(si,5) == 0
            continue
        end
        
        from_pole = abs(cell_struct(one_cell_id).Spots(si,3));
        if (.5*cell_struct(one_cell_id).Cell_Y_Axis+2-from_pole)/(.5*cell_struct(one_cell_id).Cell_Y_Axis) > -5
            
            cell_center_x = [cell_center_x cell_struct(one_cell_id).Spots(si,2)];
            cell_center_y = [cell_center_y cell_struct(one_cell_id).Spots(si,3)];
            cell_center_x2d = [cell_center_x2d spot_row];
            cell_center_y2d = [cell_center_y2d spot_col];
            cell_center_z = [cell_center_z cell_struct(one_cell_id).Spots(si,4)];
            cell_region = [cell_region cell_struct(one_cell_id).Spots(si,5)];
            if cell_struct(one_cell_id).Spots(si,5) == 4
                cell_color = [cell_color; [1 0 1]];
            elseif cell_struct(one_cell_id).Spots(si,5) == 3
                cell_color = [cell_color; [0 0 1]];
            elseif cell_struct(one_cell_id).Spots(si,5) == 1
                cell_color = [cell_color; [0 0 0]];
            elseif cell_struct(one_cell_id).Spots(si,5) == 2
                cell_color = [cell_color; [0 1 0]];
            end
            
            [cell_region,region_sort] = sort(cell_region);
            cell_center_x = cell_center_x(region_sort);
            cell_center_y = cell_center_y(region_sort);
            cell_center_x2d = cell_center_x2d(region_sort);
            cell_center_y2d = cell_center_y2d(region_sort);
            cell_center_z = cell_center_z(region_sort);
            
            
        end
    end
    
    subplot(1,2,1)
    if sum(cell_region==3) ~= 0
        gscatter(cell_center_x,cell_center_z,cell_region,'kgbm','.',10,'filled');grid on
    else
        gscatter(cell_center_x,cell_center_z,cell_region,'kgm','.',10,'filled');grid on
    end
    title(strcat(['Cell ' num2str(one_cell_id) ' XZ Spot Projection']),'FontSize',24)
    ylabel('Vertical Axis (Z)','FontSize',24)
    xlabel('Short Axis (X)','FontSize',24)
    
    subplot(1,2,2)
    
    if regular == 1
        plot(cell_struct(one_cell_id).Transformed_Boundaries{1,1}(:,2),cell_struct(one_cell_id).Transformed_Boundaries{1,1}(:,1),'k','LineWidth',3)
        hold on
        for i = 1:length(cell_struct(one_cell_id).Transformed_Dapi_Boundaries)
            plot(cell_struct(one_cell_id).Transformed_Dapi_Boundaries{i,1}(:,2),cell_struct(one_cell_id).Transformed_Dapi_Boundaries{i,1}(:,1),'m','LineWidth',3)
            hold on
        end
    elseif expanded == 1
        plot(cell_struct(one_cell_id).Transformed_Expanded_Boundaries{1,1}(:,2),cell_struct(one_cell_id).Transformed_Expanded_Boundaries{1,1}(:,1),'k','LineWidth',3)
        hold on
        for i = 1:length(cell_struct(one_cell_id).Transformed_Dapi_Boundaries)
            plot(cell_struct(one_cell_id).Transformed_Dapi_Boundaries{i,1}(:,2),cell_struct(one_cell_id).Transformed_Dapi_Boundaries{i,1}(:,1),'m','LineWidth',3)
            hold on
        end
    elseif constricted == 1
        plot(cell_struct(one_cell_id).Transformed_Constricted_Boundaries{1,1}(:,2),cell_struct(one_cell_id).Transformed_Constricted_Boundaries{1,1}(:,1),'k','LineWidth',3)
        hold on
        for i = 1:length(cell_struct(one_cell_id).Transformed_Dapi_Boundaries)
            plot(cell_struct(one_cell_id).Transformed_Dapi_Boundaries{i,1}(:,2),cell_struct(one_cell_id).Transformed_Dapi_Boundaries{i,1}(:,1),'m','LineWidth',3)
            hold on
        end
        
    end
    
    membrane_boundaries = zeros(length(cell_struct(one_cell_id).Transformed_Boundaries{1,1}),2);
    
    
    if regular == 1
        membrane_boundaries = zeros(length(cell_struct(one_cell_id).Transformed_Boundaries{1,1}),2);
        for sk = 1:length(cell_struct(one_cell_id).Transformed_Boundaries{1,1})
            point_y = cell_struct(one_cell_id).Boundaries{1,1}(sk,1) - cell_struct(one_cell_id).Center(2);
            point_x = cell_struct(one_cell_id).Boundaries{1,1}(sk,2) - cell_struct(one_cell_id).Center(1);
            point_angle = atan(point_y/point_x);
            h = sqrt((point_y)^2+(point_x)^2);
            h_prime = h - membrane_pixels;
            sign_x = sign(point_y);
            sign_y = sign(point_x);
            x_prime = abs(h_prime*sin(point_angle));
            y_prime = abs(h_prime*cos(point_angle));
            membrane_boundaries(sk,1) = sign_y*y_prime;
            membrane_boundaries(sk,2) = sign_x*x_prime;
        end
    elseif expanded == 1
        membrane_boundaries = zeros(length(cell_struct(one_cell_id).Transformed_Expanded_Boundaries{1,1}),2);
        for sk = 1:length(cell_struct(one_cell_id).Transformed_Expanded_Boundaries{1,1})
            point_y = cell_struct(one_cell_id).Expanded_Boundaries{1,1}(sk,1) - cell_struct(one_cell_id).Center(2);
            point_x = cell_struct(one_cell_id).Expanded_Boundaries{1,1}(sk,2) - cell_struct(one_cell_id).Center(1);
            point_angle = atan(point_y/point_x);
            h = sqrt((point_y)^2+(point_x)^2);
            h_prime = h - membrane_pixels;
            sign_x = sign(point_y);
            sign_y = sign(point_x);
            x_prime = abs(h_prime*sin(point_angle));
            y_prime = abs(h_prime*cos(point_angle));
            membrane_boundaries(sk,1) = sign_y*y_prime;
            membrane_boundaries(sk,2) = sign_x*x_prime;
        end
    elseif constricted == 1
        membrane_boundaries = zeros(length(cell_struct(one_cell_id).Transformed_Constricted_Boundaries{1,1}),2);
        for sk = 1:length(cell_struct(one_cell_id).Transformed_Constricted_Boundaries{1,1})
            point_y = cell_struct(one_cell_id).Constricted_Boundaries{1,1}(sk,1) - cell_struct(one_cell_id).Center(2);
            point_x = cell_struct(one_cell_id).Constricted_Boundaries{1,1}(sk,2) - cell_struct(one_cell_id).Center(1);
            point_angle = atan(point_y/point_x);
            h = sqrt((point_y)^2+(point_x)^2);
            h_prime = h - membrane_pixels;
            sign_x = sign(point_y);
            sign_y = sign(point_x);
            x_prime = abs(h_prime*sin(point_angle));
            y_prime = abs(h_prime*cos(point_angle));
            membrane_boundaries(sk,1) = sign_y*y_prime;
            membrane_boundaries(sk,2) = sign_x*x_prime;
        end
        
    end
    
    
    new_membrane_boundaries = membrane_boundaries;
    row_mem = new_membrane_boundaries(:,2);
    col_mem = new_membrane_boundaries(:,1);
    
    new_membrane_boundaries(:,1) = (col_mem*sin(cell_struct(one_cell_id).Cell_Angle)+row_mem*cos(cell_struct(one_cell_id).Cell_Angle));
    new_membrane_boundaries(:,2) = (col_mem*cos(cell_struct(one_cell_id).Cell_Angle)-row_mem*sin(cell_struct(one_cell_id).Cell_Angle));
    
    
    plot(new_membrane_boundaries(:,2),new_membrane_boundaries(:,1),'g','LineWidth',3)
    hold on
    
    pole_x = -1*(cell_struct(one_cell_id).Cell_X_Axis-1):(cell_struct(one_cell_id).Cell_X_Axis-1);
    pole_y = (.5*cell_struct(one_cell_id).Cell_Y_Axis-pole_pixels)*ones(1,length(-1*(cell_struct(one_cell_id).Cell_X_Axis-1):(cell_struct(one_cell_id).Cell_X_Axis-1)));
    neg_pole_y = -1*pole_y;
    
    plot(pole_x,pole_y,'b','LineWidth',2)
    hold on
    plot(pole_x,neg_pole_y,'b','LineWidth',2)
    hold on
    if sum(cell_region==3) ~= 0
        gscatter(cell_center_x2d,cell_center_y2d,cell_region,'kgbm','.',10,'filled');grid on
    else
        gscatter(cell_center_x2d,cell_center_y2d,cell_region,'kgm','.',10,'filled');grid on
    end
    
    set(gcf,'position',[243,868,868,667])
    prompt={'Keep Cell(1 == Yes && 0== NO) ?'};   % A box will take in the values for the X/Ytranslation
    cell_title=strcat(['Cell ' num2str(one_cell_id) ' Check']);
    answer=inputdlg(prompt,cell_title);
    discard_cell = str2num(answer{1});
    if discard_cell == 0
        discard_id = find(all_cell_ids==one_cell_id);
        cell_region_struct(discard_id) = [];
        cell_region_struct_2d(discard_id) = [];
        all_cell_ids(discard_id) = [];
    end
    
    close all
    
    
end


