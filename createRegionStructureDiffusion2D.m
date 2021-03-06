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
field12 = 'Whole_Area';
field13 = 'Pole_Area';
field14 = 'Membrane_Area';
field15 = 'Middle_Area';
field16 = 'Cytoplasm_Area';

cell_region_struct = struct(field1,[],field2, [], field3, [], field4, [],field5, [],field6, [],field7, [],field8, [],field9, [],field10, [],field11, [],field12, [],field13, [],field14, [],field15, [],field16, []);

pole_definition = .1;
middle_definition = 0; % What percentage of the way to the end of the cell should middle axis extend (should be less than 1-pole_definition)
membrane_definition = .2; %Distance to membrane (in microns) to be considered a membrane spot

membrane_pixels = membrane_definition/.130; %.130 microns/pixel

percent_off_center = 0; % How Far off center axis (in percent of cell width) to be considered a middle pixel

rid = 1;
for ri = 1:length(cell_struct)
    if size(cell_struct(ri).Spots) == [0,0];
        
        continue
    end
    cell_region_struct(rid).Cell = ri;
    cell_region_struct(rid).Total_Spots = cell_struct(ri).Num_Spots;
    cell_region_struct(rid).Length = cell_struct(ri).Cell_Y_Axis;
    pole = 0;
    membrane = 0;
    middle = 0;
    cytoplasm = 0;
    pole_diffusion = [];
    membrane_diffusion = [];
    middle_diffusion = [];
    cytoplasm_diffusion = [];
  
    
    cell_region_struct(rid).Whole_Area = sum(mask(:)==ri);
    d_pole = .5*cell_struct(ri).Cell_Y_Axis - (1-pole_definition)*.5*cell_struct(ri).Cell_Y_Axis;
    d = .5*cell_struct(ri).Cell_X_Axis - d_pole;
    c = cell_struct(ri).Cell_X_Axis*sqrt(1-(d/(.5*cell_struct(ri).Cell_X_Axis))^2);
    theta = 2*asin(c/cell_struct(ri).Cell_X_Axis); %radians
    %pole_angle = acos((cell_struct(ri).Cell_X_Axis-d_pole)/cell_struct(ri).Cell_X_Axis);
    pole_area = 2*(((.5*cell_struct(ri).Cell_X_Axis)^2)/2)*(theta-sin(theta));
    cell_region_struct(rid).Pole_Area = 2*pole_area;
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
    
    membrane_boundaries(:,1) = membrane_boundaries(:,1) + cell_struct(ri).Center(1);
    membrane_boundaries(:,2) = membrane_boundaries(:,2) + cell_struct(ri).Center(2);
    
    membrane_boundaries = int32(membrane_boundaries);
    
    for jid = 1:length(membrane_boundaries(:,1))-1
        if abs(membrane_boundaries(jid+1,1)-membrane_boundaries(jid,1)) > 1 || abs(membrane_boundaries(jid+1,2)-membrane_boundaries(jid,2)) > 1
            line = makeLine([membrane_boundaries(jid+1,2),membrane_boundaries(jid,2)],[membrane_boundaries(jid+1,1),membrane_boundaries(jid,1)]);
            membrane_boundaries = [membrane_boundaries; line];
            
        end
    end
    
    i_cell = zeros(size(mask));
    i_membrane = zeros(size(mask));
    i_not_membane = zeros(size(mask));
    
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
    
    cell_region_struct(rid).Membrane_Area = sum(i_not_membrane(:)==1);
    cell_region_struct(rid).Middle_Area = middle_definition*cell_struct(ri).Cell_Y_Axis*percent_off_center*cell_struct(ri).Cell_X_Axis;
    
    
    cell_region_struct(rid).Cytoplasm_Area = max([0,cell_region_struct(rid).Whole_Area - cell_region_struct(rid).Pole_Area - cell_region_struct(rid).Membrane_Area - cell_region_struct(rid).Middle_Area]);
    
    for sri = 1:length(cell_struct(ri).Spots(:,1))
        spot_row = int32(cell_struct(ri).Spots(sri,7));
        spot_col = int32(cell_struct(ri).Spots(sri,8));
        spot_z = int32(cell_struct(ri).Spots(sri,9));
        from_pole = abs(cell_struct(ri).Spots(sri,3));
        
        
        
        if from_pole > (1-pole_definition)*.5*cell_struct(ri).Cell_Y_Axis 
            cell_struct(ri).Spots(sri,5) = 3;
            pole = pole+1;
            pole_diffusion = [pole_diffusion spot_struct(cell_struct(ri).Spots(sri,1)).DiffusionCoefficient];
            continue
        elseif sum(membrane_row==spot_row & membrane_col==spot_col) == 1
            cell_struct(ri).Spots(sri,5) = 1;
            membrane = membrane + 1;
            membrane_diffusion = [membrane_diffusion spot_struct(cell_struct(ri).Spots(sri,1)).DiffusionCoefficient];
            continue
        elseif abs(spot_struct(cell_struct(ri).Spots(sri,1)).Collapsed_2D_Coordinate(1)) < percent_off_center*.5*cell_struct(ri).Cell_X_Axis && abs(spot_struct(cell_struct(ri).Spots(sri,1)).Collapsed_2D_Coordinate(2)) < middle_definition*.5*cell_struct(ri).Cell_Y_Axis
            cell_struct(ri).Spots(sri,5) = 4;
            middle = middle + 1;
            middle_diffusion = [middle_diffusion spot_struct(cell_struct(ri).Spots(sri,1)).DiffusionCoefficient];
            continue
        else
            cell_struct(ri).Spots(sri,5) = 2;
            cytoplasm = cytoplasm + 1;
            cytoplasm_diffusion = [cytoplasm_diffusion spot_struct(cell_struct(ri).Spots(sri,1)).DiffusionCoefficient];
            continue
        end
    end
         
    cell_region_struct(rid).Pole_Spots = pole;
    cell_region_struct(rid).Membrane_Spots = membrane;
    cell_region_struct(rid).Middle_Spots = middle;
    cell_region_struct(rid).Cytoplasm_Spots = cytoplasm;
    cell_region_struct(rid).Pole_Diffusion = pole_diffusion;
    cell_region_struct(rid).Membrane_Diffusion = membrane_diffusion;
    cell_region_struct(rid).Middle_Diffusion = middle_diffusion;
    cell_region_struct(rid).Cytoplasm_Diffusion = cytoplasm_diffusion;
    rid = rid + 1;
    
end
    

