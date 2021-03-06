mem_correct = [0,.1,.25,];
mem_def = [.15,.20,.25];
pole_def = [.1,.5,.75];

params = zeros(3,27);

mem_rich = [];
pole_rich = [];
cyto_rich = [];
mid_rich = [];
mem_spots = [];
pole_spots = [];
cyto_spots = [];
mid_spots = [];
mem_vol = [];
pole_vol = [];
cyto_vol = [];
mid_vol = [];

jid = 1;

for mc = 1:3
    for md = 1:3
        for pd = 1:3
             
            clear cell_region_struct
            field1 = 'Cell'; % All Objects, single and multi, labeled
            field2 = 'Total_Spots';
            field3 = 'Pole_Spots';
            field4 = 'Membrane_Spots';
            field5 = 'Pole_plus_Membrane_Spots';
            field6 = 'Middle_Spots';
            field7 = 'Cytoplasm_Spots';
            field8 = 'Length';
            field9 = 'Whole_Volume';
            field10 = 'Pole_Volume';
            field11 = 'Membrane_Volume';
            field12 = 'Pole_Plus_Membrane_Volume';
            field13 = 'Middle_Volume';
            field14 = 'Cytoplasm_Volume';
            
            cell_region_struct = struct(field1,[],field2, [], field3, [], field4, [],field5, [],field6, []);
            
            %Pick which type of boundaries (regular, expanded, or constricted) to use
            % based on the alignment from boundary_histogram.m
            expanded = 0;
            regular = 1;
            constricted = 0;
            
            
            pole_definition = pole_def(pd); %microns
            middle_definition = .9; % What percentage of the way to the end of the cell should middle axis extend (should be less than 1-pole_definition)
            membrane_definition = mem_def(md); %Distance to membrane (in microns) to be considered a membrane spot
            
            membrane_pixels = membrane_definition/.130; %.130 microns/pixel
            membrane_correction = mem_correct(mc); %
            pole_pixels = pole_definition/.130;
            
            height_cutoff = .4; %In microns
            
            percent_off_center = .1; % How Far off center axis (in percent of cell width) to be considered a middle pixel
            
            rid = 1;
            for ri = 1:length(cell_struct)
                if size(cell_struct(ri).Spots) == [0,0];
                    
                    continue
                end
                cell_region_struct(rid).Cell = ri;
                %cell_region_struct(rid).Total_Spots = cell_struct(ri).Num_Spots;
                cell_region_struct(rid).Length = cell_struct(ri).Cell_Y_Axis;
                pole = 0;
                membrane = 0;
                middle = 0;
                cytoplasm = 0;
                
                cell_region_struct(rid).Whole_Volume = 0;
                v  = 0;
                x_dist_array = [];
                
                if regular == 1
                    for sk = 2:length(cell_struct(ri).Transformed_Boundaries{1,1})
                        x_dist = .130*abs(mean([cell_struct(ri).Transformed_Boundaries{1,1}(sk,2),cell_struct(ri).Transformed_Boundaries{1,1}(sk-1,2)]));
                        x_dist_array = [x_dist_array x_dist];
                        if x_dist > height_cutoff
                            angle_perc = (asin(height_cutoff/x_dist)/(pi/2));
                            cell_region_struct(rid).Whole_Volume = cell_region_struct(rid).Whole_Volume + (abs(cell_struct(ri).Transformed_Boundaries{1,1}(sk,1)-cell_struct(ri).Transformed_Boundaries{1,1}(sk-1,1))-membrane_correction)*(angle_perc*(1/2)*(pi*abs(mean([cell_struct(ri).Transformed_Boundaries{1,1}(sk,2),cell_struct(ri).Transformed_Boundaries{1,1}(sk-1,2)]))^2));
                        else
                            cell_region_struct(rid).Whole_Volume = cell_region_struct(rid).Whole_Volume + (abs(cell_struct(ri).Transformed_Boundaries{1,1}(sk,1)-cell_struct(ri).Transformed_Boundaries{1,1}(sk-1,1))-membrane_correction)*((1/2)*(pi*abs(mean([cell_struct(ri).Transformed_Boundaries{1,1}(sk,2),cell_struct(ri).Transformed_Boundaries{1,1}(sk-1,2)]))^2));
                            
                        end
                    end
                elseif expanded == 1
                    for sk = 2:length(cell_struct(ri).Transformed_Expanded_Boundaries{1,1})
                        x_dist = .130*abs(mean([cell_struct(ri).Transformed_Expanded_Boundaries{1,1}(sk,2),cell_struct(ri).Transformed_Expanded_Boundaries{1,1}(sk-1,2)]));
                        x_dist_array = [x_dist_array x_dist];
                        if x_dist > height_cutoff
                            angle_perc = (asin(height_cutoff/x_dist)/(pi/2));
                            cell_region_struct(rid).Whole_Volume = cell_region_struct(rid).Whole_Volume + (abs(cell_struct(ri).Transformed_Expanded_Boundaries{1,1}(sk,1)-cell_struct(ri).Transformed_Expanded_Boundaries{1,1}(sk-1,1))-membrane_correction)*(angle_perc*(1/2)*(pi*abs(mean([cell_struct(ri).Transformed_Expanded_Boundaries{1,1}(sk,2),cell_struct(ri).Transformed_Expanded_Boundaries{1,1}(sk-1,2)]))^2));
                        else
                            cell_region_struct(rid).Whole_Volume = cell_region_struct(rid).Whole_Volume + (abs(cell_struct(ri).Transformed_Expanded_Boundaries{1,1}(sk,1)-cell_struct(ri).Transformed_Expanded_Boundaries{1,1}(sk-1,1))-membrane_correction)*((1/2)*(pi*abs(mean([cell_struct(ri).Transformed_Expanded_Boundaries{1,1}(sk,2),cell_struct(ri).Transformed_Expanded_Boundaries{1,1}(sk-1,2)]))^2));
                            
                        end
                    end
                    
                elseif constricted == 1
                    for sk = 2:length(cell_struct(ri).Transformed_Constricted_Boundaries{1,1})
                        x_dist = .130*abs(mean([cell_struct(ri).Transformed_Constricted_Boundaries{1,1}(sk,2),cell_struct(ri).Transformed_Constricted_Boundaries{1,1}(sk-1,2)]));
                        x_dist_array = [x_dist_array x_dist];
                        if x_dist > height_cutoff
                            angle_perc = (asin(height_cutoff/x_dist)/(pi/2));
                            cell_region_struct(rid).Whole_Volume = cell_region_struct(rid).Whole_Volume + (abs(cell_struct(ri).Transformed_Constricted_Boundaries{1,1}(sk,1)-cell_struct(ri).Transformed_Constricted_Boundaries{1,1}(sk-1,1))-membrane_correction)*(angle_perc*(1/2)*(pi*abs(mean([cell_struct(ri).Transformed_Constricted_Boundaries{1,1}(sk,2),cell_struct(ri).Transformed_Constricted_Boundaries{1,1}(sk-1,2)]))^2));
                        else
                            cell_region_struct(rid).Whole_Volume = cell_region_struct(rid).Whole_Volume + (abs(cell_struct(ri).Transformed_Constricted_Boundaries{1,1}(sk,1)-cell_struct(ri).Transformed_Constricted_Boundaries{1,1}(sk-1,1))-membrane_correction)*((1/2)*(pi*abs(mean([cell_struct(ri).Transformed_Constricted_Boundaries{1,1}(sk,2),cell_struct(ri).Transformed_Constricted_Boundaries{1,1}(sk-1,2)]))^2));
                            
                        end
                    end
                    
                end
                
                cell_region_struct(rid).Middle_Volume = 0;
                pole_middle_volume = 0;
                for ji = 1:length(cell_struct(ri).Transformed_Dapi_Boundaries)
                    for jid = 2:length(cell_struct(ri).Transformed_Dapi_Boundaries{ji,1})
                        y_dist = abs(mean([cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid,1),cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid-1,1)]));
                        x_dist = .130*mean([abs(cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid,2)),abs(cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid-1,2))]);
                        
                        if (y_dist > .5*cell_struct(ri).Cell_Y_Axis - pole_pixels)
                            if x_dist > height_cutoff
                                angle_perc = (asin(height_cutoff/x_dist)/(pi/2));
                                pole_middle_volume = pole_middle_volume + abs(cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid,1)-cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid-1,1))*(angle_perc*(1/2)*(pi*mean([abs(cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid,2)),abs(cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid-1,2))])^2));
                            else
                                pole_middle_volume = pole_middle_volume + abs(cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid,1)-cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid-1,1))*((1/2)*(pi*mean([abs(cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid,2)),abs(cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid-1,2))])^2));
                            end
                        end
                        
                        if x_dist > height_cutoff
                            angle_perc = (asin(height_cutoff/x_dist)/(pi/2));
                            cell_region_struct(rid).Middle_Volume = cell_region_struct(rid).Middle_Volume + abs(cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid,1)-cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid-1,1))*(angle_perc*(1/2)*(pi*mean([abs(cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid,2)),abs(cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid-1,2))])^2));
                        else
                            cell_region_struct(rid).Middle_Volume = cell_region_struct(rid).Middle_Volume + abs(cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid,1)-cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid-1,1))*((1/2)*(pi*mean([abs(cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid,2)),abs(cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid-1,2))])^2));
                        end
                    end
                end
                
                
                
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
                
                not_membrane_volume = 0;
                pole_membrane_volume = 0;
                
                check = [];
                for sk = 2:length(membrane_boundaries(:,1))
                    
                    y_dist = abs(mean([new_membrane_boundaries(sk,1),new_membrane_boundaries(sk-1,1)]));
                    x_dist = .130*abs(mean([abs(new_membrane_boundaries(sk,2)),abs(new_membrane_boundaries(sk-1,2))]));
                    
                    
                    if (y_dist > .5*cell_struct(ri).Cell_Y_Axis - pole_pixels)
                        if x_dist > height_cutoff
                            angle_perc = (asin(height_cutoff/x_dist)/(pi/2));
                            pole_membrane_volume = pole_membrane_volume + abs(new_membrane_boundaries(sk,1)-new_membrane_boundaries(sk-1,1))*(angle_perc*(1/2)*(pi*mean([abs(new_membrane_boundaries(sk,2)),abs(new_membrane_boundaries(sk-1,2))])^2));
                            
                        else
                            pole_membrane_volume = pole_membrane_volume + abs(new_membrane_boundaries(sk,1)-new_membrane_boundaries(sk-1,1))*((1/2)*(pi*mean([abs(new_membrane_boundaries(sk,2)),abs(new_membrane_boundaries(sk-1,2))])^2));
                        end
                    end
                    
                    
                    
                    if x_dist > height_cutoff
                        angle_perc = (asin(height_cutoff/x_dist)/(pi/2));
                        not_membrane_volume = not_membrane_volume + abs(new_membrane_boundaries(sk,1)-new_membrane_boundaries(sk-1,1))*(angle_perc*(1/2)*(pi*mean([abs(new_membrane_boundaries(sk,2)),abs(new_membrane_boundaries(sk-1,2))])^2));
                        
                    else
                        not_membrane_volume = not_membrane_volume + abs(new_membrane_boundaries(sk,1)-new_membrane_boundaries(sk-1,1))*((1/2)*(pi*mean([abs(new_membrane_boundaries(sk,2)),abs(new_membrane_boundaries(sk-1,2))])^2));
                    end
                    check = [check abs(new_membrane_boundaries(sk,1)-new_membrane_boundaries(sk-1,1))*((1/2)*(pi*mean([abs(new_membrane_boundaries(sk,2)),abs(new_membrane_boundaries(sk-1,2))])^2))];
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
                
                %cell_region_struct(rid).Cytoplasm_Volume = max([0,cell_region_struct(rid).Whole_Volume - cell_region_struct(rid).Pole_Volume - cell_region_struct(rid).Membrane_Volume - cell_region_struct(rid).Middle_Volume]);
                cell_region_struct(rid).Cytoplasm_Volume = max([0, not_membrane_volume - cell_region_struct(rid).Middle_Volume]);
                
                dapi_mask_i = d4 == ri;
                [dapi_row, dapi_col] = ind2sub(size(dapi_mask_i),find(dapi_mask_i));
                
                
                
                for sri = 1:length(cell_struct(ri).Spots(:,1))
                    %cell_spot_coords(qr,:) = [gs spot_struct(gs).Transform_3D_Coordinate(1) spot_struct(gs).Transform_3D_Coordinate(2) spot_struct(gs).Transform_3D_Coordinate(3) spot_struct(gs).Region spot_struct(gs).Coordinate(1) spot_struct(gs).Coordinate(2) spot_struct(gs).Coordinate(3)];
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
                    
                    if sum(membrane_row==new_spot_row & membrane_col==new_spot_col) == 1
                        cell_struct(ri).Spots(sri,5) = 1;
                        membrane = membrane + 1;
                        continue
                    elseif sum(dapi_row==new_spot_row & dapi_col==new_spot_col) == 1
                        middle = middle + 1;
                        cell_struct(ri).Spots(sri,5) = 4;
                        continue
                    elseif (from_pole <= .5*cell_struct(ri).Cell_Y_Axis - pole_pixels) && sum(cytoplasm_row==new_spot_row & cytoplasm_col==new_spot_col) == 1
                        cytoplasm = cytoplasm + 1;
                        cell_struct(ri).Spots(sri,5) = 2;
                        continue
                    elseif (from_pole > .5*cell_struct(ri).Cell_Y_Axis - pole_pixels) && (sum(cytoplasm_row==new_spot_row & cytoplasm_col==new_spot_col) == 1)
                        pole = pole+1;
                        cell_struct(ri).Spots(sri,5) = 3;
                        continue
                    else
                        cell_struct(ri).Spots(sri,5) = 0;
                    end
                end
                
                cell_region_struct(rid).Pole_Spots = pole;
                cell_region_struct(rid).Membrane_Spots = membrane;
                cell_region_struct(rid).Middle_Spots = middle;
                cell_region_struct(rid).Cytoplasm_Spots = cytoplasm;
                cell_region_struct(rid).Total_Spots = cytoplasm + pole + membrane + middle;
                rid = rid + 1;
                
            end
            
            membrane_enrichment = [];
            cytoplasm_enrichment = [];
            pole_enrichment = [];
            middle_enrichment = [];
            membrane_spots = [];
            cytoplasm_spots = [];
            pole_spots = [];
            middle_spots = [];
            membrane_volume = [];
            cytoplasm_volume = [];
            pole_volume = [];
            middle_volume = [];
            
            
            
            
            for id = 1:length(cell_region_struct)
                
                if id == 4
                    continue
                end
                
                if cell_region_struct(id).Membrane_Volume == 0
                    continue
                end
                
                if cell_region_struct(id).Cytoplasm_Volume == 0
                    continue
                end
                
                if cell_region_struct(id).Pole_Volume == 0
                    continue
                end
                
                if cell_region_struct(id).Middle_Volume == 0
                    continue
                end
                
                membrane_enrichment = [membrane_enrichment (cell_region_struct(id).Membrane_Spots/cell_region_struct(id).Total_Spots)/(cell_region_struct(id).Membrane_Volume/cell_region_struct(id).Whole_Volume)];
                
                cytoplasm_enrichment = [cytoplasm_enrichment (cell_region_struct(id).Cytoplasm_Spots/cell_region_struct(id).Total_Spots)/(cell_region_struct(id).Cytoplasm_Volume/cell_region_struct(id).Whole_Volume)];
                
                pole_enrichment = [pole_enrichment (cell_region_struct(id).Pole_Spots/cell_region_struct(id).Total_Spots)/(cell_region_struct(id).Pole_Volume/cell_region_struct(id).Whole_Volume)];
                
                middle_enrichment = [middle_enrichment (cell_region_struct(id).Middle_Spots/cell_region_struct(id).Total_Spots)/(cell_region_struct(id).Middle_Volume/cell_region_struct(id).Whole_Volume)];
                
                membrane_spots = [membrane_spots cell_region_struct(id).Membrane_Spots];
                
                cytoplasm_spots = [cytoplasm_spots cell_region_struct(id).Cytoplasm_Spots];
                
                pole_spots = [pole_spots cell_region_struct(id).Pole_Spots];
                
                middle_spots = [middle_spots cell_region_struct(id).Middle_Spots];
                
                membrane_volume = [membrane_volume cell_region_struct(id).Membrane_Volume];
                
                cytoplasm_volume = [cytoplasm_volume cell_region_struct(id).Cytoplasm_Volume];
                
                pole_volume = [pole_volume cell_region_struct(id).Pole_Volume];
                
                middle_volume = [middle_volume cell_region_struct(id).Middle_Volume];
                
            end
            
            
            enrichment = [membrane_enrichment;cytoplasm_enrichment;pole_enrichment;middle_enrichment];
            mean_enrichment = [mean(membrane_enrichment);mean(cytoplasm_enrichment);mean(pole_enrichment);mean(middle_enrichment)];

            mem_rich = [mem_rich mean(membrane_enrichment)];
            pole_rich = [pole_rich mean(pole_enrichment)];
            cyto_rich = [cyto_rich mean(cytoplasm_enrichment)];
            mid_rich = [mid_rich mean(middle_enrichment)];
            mem_spots = [mem_spots mean(membrane_spots)];
            pole_spots = [pole_spots mean(pole_spots)];
            cyto_spots = [cyto_spots mean(cytoplasm_spots)];
            mid_spots = [mid_spots mean(middle_spots)];
            mem_vol = [mem_vol mean(membrane_volume)];
            pole_vol = [pole_vol mean(pole_volume)];
            cyto_vol = [cyto_vol mean(cytoplasm_volume)];
            mid_vol = [];
            
            
            
            
        end
    end
end

jid = 1;

for mc = 1:3
    for md = 1:3
        for pd = 1:3
            params(1,jid) = mem_correct(mc);
            params(2,jid) = mem_def(md);
            params(3,jid) = pole_def(pd);
            jid = jid + 1;
        end
    end
end

mc_md = [];
mc_pd = [];
md_pd = [];

for j = 1:27
    if params(3,j) == .5
        mc_md = [mc_md j];
    end
    if params(2,j) == .25
        mc_pd = [mc_pd j];
    end
    if params(1,j) == .1
        md_pd = [md_pd j];
    end
end