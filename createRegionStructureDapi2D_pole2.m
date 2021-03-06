% clear cell_region_struct_2d
field1 = 'Cell'; % All Objects, single and multi, labeled
field2 = 'Total_Spots';
field3 = 'Pole_Spots';
field4 = 'Membrane_Spots';
field5 = 'Pole_plus_Membrane_Spots';
field6 = 'Middle_Spots';
field7 = 'Cytoplasm_Spots';
field8 = 'Length';
field9 = 'Whole_Area';
field10 = 'Pole_Area';
field11 = 'Membrane_Area';
field12 = 'Pole_Plus_Membrane_Area';
field13 = 'Middle_Area';
field14 = 'Cytoplasm_Area';

cell_region_struct_2d = struct(field1,[],field2, [], field3, [], field4, [],field5, [],field6, []);

%Pick which type of boundaries (regular, expanded, or constricted) to use
% based on the alignment from boundary_histogram.m
expanded = 0;
regular = 1;
constricted = 0;



middle_definition = .9; % What percentage of the way to the end of the cell should middle axis extend (should be less than 1-pole_definition)
membrane_definition = (boundary/2)*.130; %Distance to membrane (in microns) to be considered a membrane spot
pole_definition = membrane_definition+.1; %microns

membrane_pixels = membrane_definition/.130; %.130 microns/pixel
membrane_correction = .1*boundary; %
pole_pixels = pole_definition/.130;

height_cutoff = .65; %In microns
percent_off_center = .1; % How Far off center axis (in percent of cell width) to be considered a middle pixel
middle_section_cutoff = 100; %In nm, how big is the middle section to be considered 2D? 100 = 50 nm above and below 0

up_check = [];
down_check = [];

rid = 1;
for ri = 1:length(cell_struct)
    if size(cell_struct(ri).Spots) == [0,0] | isempty(cell_struct(ri).Dapi_Boundaries);
        
        continue
    end
    cell_region_struct_2d(rid).Cell = ri;
    
    cell_region_struct_2d(rid).Length = cell_struct(ri).Cell_Y_Axis;
    pole = 0;
    membrane = 0;
    middle = 0;
    cytoplasm = 0;
    
    cell_region_struct_2d(rid).Whole_Area = 0;
    v  = 0;
    x_dist_array = [];
    
    corrected_boundaries = zeros(length(cell_struct(ri).Boundaries{1,1}),2);
    for sk = 1:length(cell_struct(ri).Transformed_Boundaries{1,1})
        point_y = cell_struct(ri).Boundaries{1,1}(sk,1) - cell_struct(ri).Center(2);
        point_x = cell_struct(ri).Boundaries{1,1}(sk,2) - cell_struct(ri).Center(1);
        point_angle = atan(point_y/point_x);
        h = sqrt((point_y)^2+(point_x)^2);
        h_prime = h - membrane_correction;
        sign_x = sign(point_y);
        sign_y = sign(point_x);
        x_prime = abs(h_prime*sin(point_angle));
        y_prime = abs(h_prime*cos(point_angle));
        corrected_boundaries(sk,1) = sign_y*y_prime;
        corrected_boundaries(sk,2) = sign_x*x_prime;
    end
    
    if regular == 1
        %cell_region_struct_2d(rid).Whole_Area = polyarea(cell_struct(ri).Boundaries{1,1}(:,1),cell_struct(ri).Boundaries{1,1}(:,2));
        cell_region_struct_2d(rid).Whole_Area = polyarea(corrected_boundaries(:,1),corrected_boundaries(:,2));
    elseif expanded == 1
        cell_region_struct_2d(rid).Whole_Area = polyarea(cell_struct(ri).Expanded_Boundaries{1,1}(:,1),cell_struct(ri).Expanded_Boundaries{1,1}(:,2));
        
    elseif constricted == 1
        cell_region_struct_2d(rid).Whole_Area = polyarea(cell_struct(ri).Constricted_Boundaries{1,1}(:,1),cell_struct(ri).Constricted_Boundaries{1,1}(:,2));
        
    end
    
    cell_region_struct_2d(rid).Middle_Area = 0;
    pole_middle_volume = 0;
    dapi_x_shift = [];
    for ji = 1:length(cell_struct(ri).Transformed_Dapi_Boundaries)
        
        for jid = 2:length(cell_struct(ri).Transformed_Dapi_Boundaries{ji,1})
            dapi_x_shift = [dapi_x_shift cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid,2)];
        end
    end
    
    dapi_x_shift = mean(dapi_x_shift);
    middle_area = 0;
    %dapi_x_shift = 0;
    for ji = 1:length(cell_struct(ri).Transformed_Dapi_Boundaries)
        cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(:,2) = cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(:,2) - dapi_x_shift;
        middle_area = middle_area + polyarea(cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(:,1),cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(:,2));
    end
    
    
    cell_region_struct_2d(rid).Middle_Area = middle_area;
    

    if regular == 1
        membrane_boundaries = zeros(length(cell_struct(ri).Transformed_Boundaries{1,1}),2);
        for sk = 1:length(cell_struct(ri).Transformed_Boundaries{1,1})
            point_y = cell_struct(ri).Boundaries{1,1}(sk,1) - cell_struct(ri).Center(2);
            point_x = cell_struct(ri).Boundaries{1,1}(sk,2) - cell_struct(ri).Center(1);
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
        membrane_boundaries = zeros(length(cell_struct(ri).Transformed_Expanded_Boundaries{1,1}),2);
        for sk = 1:length(cell_struct(ri).Transformed_Expanded_Boundaries{1,1})
            point_y = cell_struct(ri).Expanded_Boundaries{1,1}(sk,1) - cell_struct(ri).Center(2);
            point_x = cell_struct(ri).Expanded_Boundaries{1,1}(sk,2) - cell_struct(ri).Center(1);
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
        membrane_boundaries = zeros(length(cell_struct(ri).Transformed_Constricted_Boundaries{1,1}),2);
        for sk = 1:length(cell_struct(ri).Transformed_Constricted_Boundaries{1,1})
            point_y = cell_struct(ri).Constricted_Boundaries{1,1}(sk,1) - cell_struct(ri).Center(2);
            point_x = cell_struct(ri).Constricted_Boundaries{1,1}(sk,2) - cell_struct(ri).Center(1);
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
    
    new_membrane_boundaries(:,1) = (col_mem*sin(cell_struct(ri).Cell_Angle)+row_mem*cos(cell_struct(ri).Cell_Angle));
    new_membrane_boundaries(:,2) = (col_mem*cos(cell_struct(ri).Cell_Angle)-row_mem*sin(cell_struct(ri).Cell_Angle));
    
    
    pole_middle_up = [];
    pole_middle_down = [];
    pole_membrane_up = [];
    pole_membrane_down = [];
    
    pole_membrane_up = new_membrane_boundaries((new_membrane_boundaries(:,1)>.5*cell_struct(ri).Cell_Y_Axis - pole_pixels),:);
    if (new_membrane_boundaries(1,1)>0.5*cell_struct(ri).Cell_Y_Axis - pole_pixels)
        mem_switch = 1;
    else
        mem_switch = -1;
    end
    pre_mem_switch = mem_switch;
    
    pmu_indices = [];
    
    for pm = 1:length(new_membrane_boundaries)
        pre_mem_switch = mem_switch;
        if (new_membrane_boundaries(pm,1)>0.5*cell_struct(ri).Cell_Y_Axis - pole_pixels)
            mem_switch = 1;
        else
            mem_switch = -1;
        end
        
        if mem_switch ~= pre_mem_switch
            pmu_indices = [pmu_indices; pm-1; pm];
        end
    end
    
    up_ind=find((new_membrane_boundaries(:,1)>0.5*cell_struct(ri).Cell_Y_Axis - pole_pixels));
    pole_membrane_up_indices_min = pmu_indices(1);
    pole_membrane_up_indices_max = pmu_indices(4);
    pole_membrane_up_indices_min_one = pmu_indices(2);
    pole_membrane_up_indices_max_one = pmu_indices(3);
    
    up_slope_min = (new_membrane_boundaries(pole_membrane_up_indices_min,1)-new_membrane_boundaries(pole_membrane_up_indices_min_one,1))/(new_membrane_boundaries(pole_membrane_up_indices_min,2)-new_membrane_boundaries(pole_membrane_up_indices_min_one,2));
    up_slope_max = (new_membrane_boundaries(pole_membrane_up_indices_max,1)-new_membrane_boundaries(pole_membrane_up_indices_max_one,1))/(new_membrane_boundaries(pole_membrane_up_indices_max,2)-new_membrane_boundaries(pole_membrane_up_indices_max_one,2));
    intercept_up_min = new_membrane_boundaries(pole_membrane_up_indices_min,1) - (up_slope_min*new_membrane_boundaries(pole_membrane_up_indices_min,2));
    intercept_up_max = new_membrane_boundaries(pole_membrane_up_indices_max,1) - (up_slope_max*new_membrane_boundaries(pole_membrane_up_indices_max,2));
    up_pole_x_min = [(.5*cell_struct(ri).Cell_Y_Axis - pole_pixels), ((.5*cell_struct(ri).Cell_Y_Axis - pole_pixels)-intercept_up_min)/up_slope_min];
    up_pole_x_max = [(.5*cell_struct(ri).Cell_Y_Axis - pole_pixels), ((.5*cell_struct(ri).Cell_Y_Axis - pole_pixels)-intercept_up_max)/up_slope_max];
    if pmu_indices(1) < min(up_ind)
        new_up_indices = [min(up_ind):max(up_ind)];
    else
        new_up_indices = [pmu_indices(1):-1:min(up_ind),max(up_ind):-1:pmu_indices(4)];
    end
       
    pole_membrane_up = [up_pole_x_min;new_membrane_boundaries(new_up_indices,:);up_pole_x_max];
    
    
    pole_membrane_down = new_membrane_boundaries((new_membrane_boundaries(:,1)<(-0.5)*cell_struct(ri).Cell_Y_Axis + pole_pixels),:);
    
    if (new_membrane_boundaries(1,1)<(-0.5)*cell_struct(ri).Cell_Y_Axis + pole_pixels)
        mem_switch = 1;
    else
        mem_switch = -1;
    end
    pre_mem_switch = mem_switch;
    
    pmd_indices = [];
    
    for pm = 1:length(new_membrane_boundaries)
        pre_mem_switch = mem_switch;
        if (new_membrane_boundaries(pm,1)<(-0.5)*cell_struct(ri).Cell_Y_Axis + pole_pixels)
            mem_switch = 1;
        else
            mem_switch = -1;
        end
        
        if mem_switch ~= pre_mem_switch
            pmd_indices = [pmd_indices; pm-1; pm];
        end
    end
    
    down_ind=find((new_membrane_boundaries(:,1)<(-0.5)*cell_struct(ri).Cell_Y_Axis + pole_pixels));
    pole_membrane_down_indices_min = pmd_indices(1);
    pole_membrane_down_indices_max = pmd_indices(4);
    pole_membrane_down_indices_min_one = pmd_indices(2);
    pole_membrane_down_indices_max_one = pmd_indices(3);
    down_slope_min = (new_membrane_boundaries(pole_membrane_down_indices_min,1)-new_membrane_boundaries(pole_membrane_down_indices_min_one,1))/(new_membrane_boundaries(pole_membrane_down_indices_min,2)-new_membrane_boundaries(pole_membrane_down_indices_min_one,2));
    down_slope_max = (new_membrane_boundaries(pole_membrane_down_indices_max,1)-new_membrane_boundaries(pole_membrane_down_indices_max_one,1))/(new_membrane_boundaries(pole_membrane_down_indices_max,2)-new_membrane_boundaries(pole_membrane_down_indices_max_one,2));
    intercept_down_min = new_membrane_boundaries(pole_membrane_down_indices_min,1) - (down_slope_min*new_membrane_boundaries(pole_membrane_down_indices_min,2));
    intercept_down_max = new_membrane_boundaries(pole_membrane_down_indices_max,1) - (down_slope_max*new_membrane_boundaries(pole_membrane_down_indices_max,2));
    down_pole_x_min = [((-0.5)*cell_struct(ri).Cell_Y_Axis + pole_pixels), (((-0.5)*cell_struct(ri).Cell_Y_Axis + pole_pixels)-intercept_down_min)/down_slope_min];
    down_pole_x_max = [((-0.5)*cell_struct(ri).Cell_Y_Axis + pole_pixels), (((-0.5)*cell_struct(ri).Cell_Y_Axis + pole_pixels)-intercept_down_max)/down_slope_max];
    if pmd_indices(1) < min(down_ind)
        new_down_indices = [min(down_ind):max(down_ind)];
    else
        new_down_indices = [pmd_indices(1):1:min(down_ind),max(down_ind):-1:pmd_indices(4)];
    end
    pole_membrane_down = [down_pole_x_min;new_membrane_boundaries(new_down_indices,:);down_pole_x_max];
    
    
    for ji = 1:length(cell_struct(ri).Transformed_Dapi_Boundaries)
        pole_middle_up = [pole_middle_up; cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(((cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(:,2)>.5*cell_struct(ri).Cell_Y_Axis - pole_pixels) & (cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(:,2)<.5*cell_struct(ri).Cell_Y_Axis - membrane_pixels)),:)];
        pole_middle_down = [pole_middle_down cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(((cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(:,2)<(-0.5)*cell_struct(ri).Cell_Y_Axis + pole_pixels)& (cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(:,2)>(-0.5)*cell_struct(ri).Cell_Y_Axis + membrane_pixels)),:)];
    end
    
    
    
    pole_area_up = max([0,polyarea(pole_membrane_up(:,2),pole_membrane_up(:,1)) - polyarea(pole_middle_up(:,2),pole_middle_up(:,1))]);
    pole_area_down = max([0,polyarea(pole_membrane_down(:,2),pole_membrane_down(:,1)) - polyarea(pole_middle_down(:,2),pole_middle_down(:,1))]);
    
    if pole_area_up == 0
        up_check = [up_check ri];
    end
    
    if pole_area_down == 0
        down_check = [down_check ri];
    end

    cell_region_struct_2d(rid).Pole_Area = pole_area_up+pole_area_down;
    
    cell_region_struct_2d(rid).Membrane_Area = max([0,cell_region_struct_2d(rid).Whole_Area - polyarea(new_membrane_boundaries(:,2),new_membrane_boundaries(:,1))]);
    
    cell_region_struct_2d(rid).Cytoplasm_Area = max([0,cell_region_struct_2d(rid).Whole_Area - cell_region_struct_2d(rid).Pole_Area - cell_region_struct_2d(rid).Membrane_Area - cell_region_struct_2d(rid).Middle_Area]);
    
    for sri = 1:length(cell_struct(ri).Spots(:,1))
        %cell_spot_coords(qr,:) = [gs spot_struct(gs).Transform_3D_Coordinate(1) spot_struct(gs).Transform_3D_Coordinate(2) spot_struct(gs).Transform_3D_Coordinate(3) spot_struct(gs).Region spot_struct(gs).Coordinate(1) spot_struct(gs).Coordinate(2) spot_struct(gs).Coordinate(3)];
        middle_check = 0;
        spot_id = cell_struct(ri).Spots(sri,1);
        spot_row = spot_struct(spot_id).Transform_3D_Coordinate(1);
        spot_col = spot_struct(spot_id).Transform_3D_Coordinate(2);
        spot_z = 0;%130*abs(cell_struct(ri).Spots(sri,4));
        
        p = rand;
        
        if p > 0.5
            continue
        end
          
        
        
        if spot_z > middle_section_cutoff/2
            cell_struct(ri).Spots(sri,5) = 0;
            continue
        end
        
        
        new_spot_row = int32((spot_row*sin(-1*cell_struct(ri).Cell_Angle)+spot_col*cos(-1*cell_struct(ri).Cell_Angle)) + cell_struct(ri).Center(2));
        new_spot_col = int32((spot_row*cos(-1*cell_struct(ri).Cell_Angle)-spot_col*sin(-1*cell_struct(ri).Cell_Angle)) + cell_struct(ri).Center(1));
        
       
        
        from_pole = abs(spot_col);
        
        cell_struct(ri).Spots(sri,5) = 0;
        
        if (from_pole > .5*cell_struct(ri).Cell_Y_Axis - pole_pixels) && inpolygon(spot_row,spot_col,new_membrane_boundaries(:,2),new_membrane_boundaries(:,1)) && middle_check == 0
            pole = pole+1;
            cell_struct(ri).Spots(sri,5) = 3;
            continue
        end
        
        if inpolygon(spot_row,spot_col,cell_struct(ri).Transformed_Boundaries{1,1}(:,2),cell_struct(ri).Transformed_Boundaries{1,1}(:,1)) && ~(inpolygon(spot_row,spot_col,new_membrane_boundaries(:,2),new_membrane_boundaries(:,1))) && middle_check == 0%sum(membrane_row==new_spot_row & membrane_col==new_spot_col) == 1
            cell_struct(ri).Spots(sri,5) = 1;
            membrane = membrane + 1;
            continue
        end
        
        for rjd = 1:length(cell_struct(ri).Transformed_Dapi_Boundaries)
            if inpolygon(spot_row,spot_col,cell_struct(ri).Transformed_Dapi_Boundaries{rjd,1}(:,2),cell_struct(ri).Transformed_Dapi_Boundaries{rjd,1}(:,1))% && inpolygon(spot_row,spot_col,new_membrane_boundaries(:,2),new_membrane_boundaries(:,1))
                middle = middle + 1;
                cell_struct(ri).Spots(sri,5) = 4;
                middle_check = 1;
                continue
            end
        end
        
        

        if (from_pole <= .5*cell_struct(ri).Cell_Y_Axis - pole_pixels) && inpolygon(spot_row,spot_col,new_membrane_boundaries(:,2),new_membrane_boundaries(:,1)) && middle_check == 0
            cytoplasm = cytoplasm + 1;
            cell_struct(ri).Spots(sri,5) = 2;
            continue
        
        else
            if middle_check == 0
                cell_struct(ri).Spots(sri,5) = 0;
            end
        end
    end
         
    cell_region_struct_2d(rid).Pole_Spots = pole;
    cell_region_struct_2d(rid).Membrane_Spots = membrane;
    cell_region_struct_2d(rid).Middle_Spots = middle;
    cell_region_struct_2d(rid).Cytoplasm_Spots = cytoplasm;
    cell_region_struct_2d(rid).Total_Spots = cytoplasm + pole + membrane + middle;
    rid = rid + 1;
    
end
    