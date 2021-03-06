% Creates a new structure called "cell_region_struct" found in workspace

clear cell_region_struct
field1 = 'Cell'; % All Objects, single and multi, labeled
field2 = 'Total_Spots';
field3 = 'Pole_Spots';
field4 = 'Pole_Diffusion';
field5 = 'Membrane_Spots';
field6 = 'Membrane_Diffusion';
field7 = 'Middle_Spots';
field8 = 'Middle_Diffusion';
field9 = 'Cytoplasm_Spots';
field10 = 'Cytoplasm_Diffusion';
field11 = 'Length';
field12 = 'Whole_Volume';
field13 = 'Pole_Volume';
field14 = 'Membrane_Volume';
field15 = 'Middle_Volume';
field16 = 'Cytoplasm_Volume';

d5 = bwlabel(d4,4);
num_dapi = max(d5(:));

for j = 1:num_dapi
    temp_dapi_mask = d5==j;
    dapi_ids = [];
    dapi_ids = unique(temp_dapi_mask.*new_mask);
    
    if sum(dapi_ids ~= 0) == 0
        d4(d5 == j) = 0;
    end
        
end

for i = 1:max(new_mask(:))
    strcat(['Working on Cell ' num2str(i)])
    
    if sum(new_mask(:)==i) == 0
        continue
    end
   
    cell_struct(i).Dapi_Boundaries = bwboundaries(d4==i);
    cell_struct(i).Transformed_Dapi_Boundaries = cell_struct(i).Dapi_Boundaries;
    
    
    if ~isempty(cell_struct(i).Dapi_Boundaries) 
        
        for di = 1:length(cell_struct(i).Dapi_Boundaries)
            
            cell_struct(i).Transformed_Dapi_Boundaries{di,1}(:,1) = cell_struct(i).Transformed_Dapi_Boundaries{di,1}(:,1) - cell_struct(i).Center(2);
            cell_struct(i).Transformed_Dapi_Boundaries{di,1}(:,2) = cell_struct(i).Transformed_Dapi_Boundaries{di,1}(:,2) - cell_struct(i).Center(1);
            dapi_row_border = cell_struct(i).Transformed_Dapi_Boundaries{di,1}(:,1);
            dapi_col_border = cell_struct(i).Transformed_Dapi_Boundaries{di,1}(:,2);
            cell_struct(i).Transformed_Dapi_Boundaries{di,1}(:,1) = (dapi_col_border*sin(cell_struct(i).Cell_Angle)+dapi_row_border*cos(cell_struct(i).Cell_Angle));
            cell_struct(i).Transformed_Dapi_Boundaries{di,1}(:,2) = (dapi_col_border*cos(cell_struct(i).Cell_Angle)-dapi_row_border*sin(cell_struct(i).Cell_Angle));
            
        end
    end
    
    
end

cell_region_struct = struct(field1,[],field2, [], field3, [], field4, [],field5, [],field6, [],field7, [],field8, [],field9, [],field10, [],field11, [],field12, [],field13, [],field14, [],field15, [],field16, []);

pole_definition = .1;
middle_definition = .9; % What percentage of the way to the end of the cell should middle axis extend (should be less than 1-pole_definition)
membrane_definition = .1; %Distance to membrane (in microns) to be considered a membrane spot

membrane_pixels = membrane_definition/.130; %.130 microns/pixel
membrane_correction = .25; %
pole_pixels = pole_definition/.130;

height_cutoff = .65; %In microns
percent_off_center = .1; % How Far off center axis (in percent of cell width) to be considered a middle pixel



rid = 1;
for ri = 1:length(cell_struct)
    if size(cell_struct(ri).Spots) == [0,0];
        
        continue
    end
    cell_region_struct(rid).Cell = ri;
    
    cell_region_struct(rid).Length = cell_struct(ri).Cell_Y_Axis;
    pole = 0;
    membrane = 0;
    middle = 0;
    cytoplasm = 0;
    pole_diffusion = [];
    membrane_diffusion = [];
    middle_diffusion = [];
    cytoplasm_diffusion = [];
    
    cell_region_struct(rid).Whole_Volume = 0;
    v  = 0;
    x_dist_array = [];
    
    for sk = 2:length(cell_struct(ri).Transformed_Boundaries{1,1})
        x_dist = .130*abs(mean([cell_struct(ri).Transformed_Boundaries{1,1}(sk,2),cell_struct(ri).Transformed_Boundaries{1,1}(sk-1,2)]));
        x_dist_array = [x_dist_array x_dist];
        if x_dist > height_cutoff
            angle_perc = (asin(height_cutoff/x_dist)/(pi/2));
            cell_region_struct(rid).Whole_Volume = cell_region_struct(rid).Whole_Volume + (abs(cell_struct(ri).Transformed_Boundaries{1,1}(sk,1)-cell_struct(ri).Transformed_Boundaries{1,1}(sk-1,1)))*(angle_perc*(1/2)*(pi*(abs(mean([cell_struct(ri).Transformed_Boundaries{1,1}(sk,2),cell_struct(ri).Transformed_Boundaries{1,1}(sk-1,2)]))-membrane_correction)^2));
        else
            cell_region_struct(rid).Whole_Volume = cell_region_struct(rid).Whole_Volume + (abs(cell_struct(ri).Transformed_Boundaries{1,1}(sk,1)-cell_struct(ri).Transformed_Boundaries{1,1}(sk-1,1)))*((1/2)*(pi*(abs(mean([cell_struct(ri).Transformed_Boundaries{1,1}(sk,2),cell_struct(ri).Transformed_Boundaries{1,1}(sk-1,2)]))-membrane_correction)^2));
            
        end
    end
    
    
    cell_region_struct(rid).Middle_Volume = 0;
    pole_middle_volume = 0;
    dapi_x_shift = [];
    for ji = 1:length(cell_struct(ri).Transformed_Dapi_Boundaries)
        
        for jid = 2:length(cell_struct(ri).Transformed_Dapi_Boundaries{ji,1})
            dapi_x_shift = [dapi_x_shift cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid,2)];
        end
    end
    
    dapi_x_shift = mean(dapi_x_shift);
    %dapi_x_shift = 0;
    for ji = 1:length(cell_struct(ri).Transformed_Dapi_Boundaries)
        cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(:,2) = cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(:,2) - dapi_x_shift;
    end
    
    for ji = 1:length(cell_struct(ri).Transformed_Dapi_Boundaries)
        for jid = 2:length(cell_struct(ri).Transformed_Dapi_Boundaries{ji,1})
            
            y_dist = abs(mean([cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid,1),cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid-1,1)]));
            x_dist = .130*(abs(mean([cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid,2),cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid-1,2)])));
            
            if (y_dist > .5*cell_struct(ri).Cell_Y_Axis - pole_pixels)
                if x_dist > height_cutoff
                    angle_perc = (asin(height_cutoff/x_dist)/(pi/2));
                    pole_middle_volume = pole_middle_volume + abs(cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid,1)-cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid-1,1))*(angle_perc*(1/2)*(pi*(abs(mean([cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid,2),cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid-1,2)])))^2));
                else
                    pole_middle_volume = pole_middle_volume + abs(cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid,1)-cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid-1,1))*((1/2)*(pi*(abs(mean([cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid,2),cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid-1,2)])))^2));
                end
            end
            
            if x_dist > height_cutoff
                angle_perc = (asin(height_cutoff/x_dist)/(pi/2));
                cell_region_struct(rid).Middle_Volume = cell_region_struct(rid).Middle_Volume + abs(cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid,1)-cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid-1,1))*(angle_perc*(1/2)*(pi*(abs(mean([cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid,2),cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid-1,2)])))^2));
            else
                cell_region_struct(rid).Middle_Volume = cell_region_struct(rid).Middle_Volume + abs(cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid,1)-cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid-1,1))*((1/2)*(pi*(abs(mean([cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid,2),cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid-1,2)])))^2));
            end
        end
    end
    
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
    
    new_membrane_boundaries = membrane_boundaries;
    row_mem = new_membrane_boundaries(:,2);
    col_mem = new_membrane_boundaries(:,1);
    
    new_membrane_boundaries(:,1) = (col_mem*sin(cell_struct(ri).Cell_Angle)+row_mem*cos(cell_struct(ri).Cell_Angle));
    new_membrane_boundaries(:,2) = (col_mem*cos(cell_struct(ri).Cell_Angle)-row_mem*sin(cell_struct(ri).Cell_Angle));
    
    not_membrane_volume = 0;
    pole_membrane_volume = 0;
    
    check = [];
    for sk = 2:length(membrane_boundaries(:,1))
        
        y_dist = abs(mean([new_membrane_boundaries(sk,1),new_membrane_boundaries(sk-1,1)]));
        x_dist = .130*abs(mean([new_membrane_boundaries(sk,2),new_membrane_boundaries(sk-1,2)]));
        
        
        if (y_dist > .5*cell_struct(ri).Cell_Y_Axis - pole_pixels)
            if x_dist > height_cutoff
                angle_perc = (asin(height_cutoff/x_dist)/(pi/2));
                pole_membrane_volume = pole_membrane_volume + abs(new_membrane_boundaries(sk,1)-new_membrane_boundaries(sk-1,1))*(angle_perc*(1/2)*(pi*abs(mean([new_membrane_boundaries(sk,2),new_membrane_boundaries(sk-1,2)]))^2));
                
            else
                pole_membrane_volume = pole_membrane_volume + abs(new_membrane_boundaries(sk,1)-new_membrane_boundaries(sk-1,1))*((1/2)*(pi*abs(mean([new_membrane_boundaries(sk,2),new_membrane_boundaries(sk-1,2)]))^2));
            end
        end
            
        
        
        if x_dist > height_cutoff
            angle_perc = (asin(height_cutoff/x_dist)/(pi/2));
            not_membrane_volume = not_membrane_volume + abs(new_membrane_boundaries(sk,1)-new_membrane_boundaries(sk-1,1))*(angle_perc*(1/2)*(pi*abs(mean([new_membrane_boundaries(sk,2),new_membrane_boundaries(sk-1,2)]))^2));
                
        else
            not_membrane_volume = not_membrane_volume + abs(new_membrane_boundaries(sk,1)-new_membrane_boundaries(sk-1,1))*((1/2)*(pi*abs(mean([new_membrane_boundaries(sk,2),new_membrane_boundaries(sk-1,2)]))^2));
        end
        check = [check abs(new_membrane_boundaries(sk,1)-new_membrane_boundaries(sk-1,1))*((1/2)*(pi*abs(mean([new_membrane_boundaries(sk,2),new_membrane_boundaries(sk-1,2)]))^2))];
    end
    
    cell_region_struct(rid).Pole_Volume = max([0,pole_membrane_volume - pole_middle_volume]);
    
    cell_region_struct(rid).Membrane_Volume = max([0,cell_region_struct(rid).Whole_Volume - not_membrane_volume]);
    
    membrane_boundaries(:,1) = membrane_boundaries(:,1) + cell_struct(ri).Center(1); %X
    membrane_boundaries(:,2) = membrane_boundaries(:,2) + cell_struct(ri).Center(2); %Y
    
    membrane_boundaries = int32(membrane_boundaries);
    
    for jid = 1:length(membrane_boundaries(:,1))-1   
        if abs(membrane_boundaries(jid+1,1)-membrane_boundaries(jid,1)) > 1 || abs(membrane_boundaries(jid+1,2)-membrane_boundaries(jid,2)) > 1
            line = makeLine([membrane_boundaries(jid+1,2),membrane_boundaries(jid,2)],[membrane_boundaries(jid+1,1),membrane_boundaries(jid,1)]);
            membrane_boundaries = [membrane_boundaries; line];
        end
    end
    
    
    
    i_cell = zeros(size(mask));
    i_membrane = zeros(size(mask));
    i_not_membrane = zeros(size(mask));
    
    for tk = 1:length(membrane_boundaries(:,1))
        
        i_membrane(membrane_boundaries(tk,2),membrane_boundaries(tk,1)) = 1;
    end
    
    for tj = 1:length(cell_struct(ri).Boundaries{1,1})
        i_cell(cell_struct(ri).Boundaries{1,1}(tj,1),cell_struct(ri).Boundaries{1,1}(tj,2)) = 1;
    end
   
    i_cell = imfill(i_cell,'holes');
    i_membrane = imfill(i_membrane,'holes');
    i_not_membrane = i_cell - i_membrane;
    [membrane_row, membrane_col] = ind2sub(size(i_not_membrane),find(i_not_membrane));     
    [cytoplasm_row, cytoplasm_col] = ind2sub(size(i_membrane),find(i_membrane));
    
    cell_region_struct(rid).Cytoplasm_Volume = max([0,cell_region_struct(rid).Whole_Volume - cell_region_struct(rid).Pole_Volume - cell_region_struct(rid).Membrane_Volume - cell_region_struct(rid).Middle_Volume]);
    %cell_region_struct(rid).Cytoplasm_Volume = max([0, not_membrane_volume - cell_region_struct(rid).Middle_Volume]);
    
    
    dapi_mask_i = d4 == ri;
    [dapi_row, dapi_col] = ind2sub(size(dapi_mask_i),find(dapi_mask_i));
    
    
    
    for sri = 1:length(cell_struct(ri).Spots(:,1))
        middle_check = 0;
        spot_id = cell_struct(ri).Spots(sri,1);
        spot_row = spot_struct(spot_id).Collapsed_2D_Coordinate(1);
        spot_col = spot_struct(spot_id).Collapsed_2D_Coordinate(2);
        spot_z = .130*abs(cell_struct(ri).Spots(sri,4));
        
        
        if spot_z > height_cutoff
            cell_struct(ri).Spots(sri,5) = 0;
            continue
        end
        
        
        new_spot_row = int32((spot_row*sin(-1*cell_struct(ri).Cell_Angle)+spot_col*cos(-1*cell_struct(ri).Cell_Angle)) + cell_struct(ri).Center(2));
        new_spot_col = int32((spot_row*cos(-1*cell_struct(ri).Cell_Angle)-spot_col*sin(-1*cell_struct(ri).Cell_Angle)) + cell_struct(ri).Center(1));
        
       
        
        from_pole = abs(spot_col);
        
        cell_struct(ri).Spots(sri,5) = 0;
        
        for rjd = 1:length(cell_struct(ri).Transformed_Dapi_Boundaries)
            if inpolygon(spot_row,spot_col,cell_struct(ri).Transformed_Dapi_Boundaries{rjd,1}(:,2),cell_struct(ri).Transformed_Dapi_Boundaries{rjd,1}(:,1))% && inpolygon(spot_row,spot_col,new_membrane_boundaries(:,2),new_membrane_boundaries(:,1))
                middle = middle + 1;
                cell_struct(ri).Spots(sri,5) = 4;
                middle_diffusion = [middle_diffusion spot_struct(cell_struct(ri).Spots(sri,1)).DiffusionCoefficient];
                middle_check = 1;
                continue
            end
        end
        
        if inpolygon(spot_row,spot_col,cell_struct(ri).Transformed_Boundaries{1,1}(:,2),cell_struct(ri).Transformed_Boundaries{1,1}(:,1)) && ~(inpolygon(spot_row,spot_col,new_membrane_boundaries(:,2),new_membrane_boundaries(:,1))) && middle_check == 0%sum(membrane_row==new_spot_row & membrane_col==new_spot_col) == 1
            cell_struct(ri).Spots(sri,5) = 1;
            membrane = membrane + 1;
            membrane_diffusion = [membrane_diffusion spot_struct(cell_struct(ri).Spots(sri,1)).DiffusionCoefficient];
            continue
            
        elseif (from_pole <= .5*cell_struct(ri).Cell_Y_Axis - pole_pixels) && inpolygon(spot_row,spot_col,new_membrane_boundaries(:,2),new_membrane_boundaries(:,1)) && middle_check == 0
            cytoplasm = cytoplasm + 1;
            cell_struct(ri).Spots(sri,5) = 2;
            cytoplasm_diffusion = [cytoplasm_diffusion spot_struct(cell_struct(ri).Spots(sri,1)).DiffusionCoefficient];
            continue
        elseif (from_pole > .5*cell_struct(ri).Cell_Y_Axis - pole_pixels) && inpolygon(spot_row,spot_col,new_membrane_boundaries(:,2),new_membrane_boundaries(:,1)) && middle_check == 0
            pole = pole+1;
            cell_struct(ri).Spots(sri,5) = 3;
            pole_diffusion = [pole_diffusion spot_struct(cell_struct(ri).Spots(sri,1)).DiffusionCoefficient];
            continue
        else
            if middle_check == 0
                cell_struct(ri).Spots(sri,5) = 0;
            end
        end
    
    end
         
    cell_region_struct(rid).Pole_Spots = pole;
    cell_region_struct(rid).Membrane_Spots = membrane;
    cell_region_struct(rid).Middle_Spots = middle;
    cell_region_struct(rid).Cytoplasm_Spots = cytoplasm;
    cell_region_struct(rid).Total_Spots = cytoplasm + pole + membrane + middle;
    
    cell_region_struct(rid).Pole_Diffusion = pole_diffusion;
    cell_region_struct(rid).Membrane_Diffusion = membrane_diffusion;
    cell_region_struct(rid).Middle_Diffusion = middle_diffusion;
    cell_region_struct(rid).Cytoplasm_Diffusion = cytoplasm_diffusion;
    rid = rid + 1;
    
end
    

