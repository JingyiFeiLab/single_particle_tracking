% Creates a new structure called "cell_region_struct" found in workspace

clear cell_region_struct
field1 = 'Cell'; % All Objects, single and multi, labeled
field2 = 'Total_Spots';
field3 = 'Membrane_Spots';
field4 = 'Middle_Spots';
field5 = 'Cytoplasm_Spots';
field6 = 'Length';
field7 = 'Whole_Volume';
field8 = 'Membrane_Volume';
field9 = 'Middle_Volume';
field10 = 'Cytoplasm_Volume';
field11 = 'Pole_Volume';

cell_region_struct = struct(field1,[],field2, [], field3, [], field4, [],field5, [],field6, [], field7, [], field8, [],field9, [],field10, []);

pole_definition = .01;
middle_definition = .9; % What percentage of the way to the end of the cell should middle axis extend (should be less than 1-pole_definition)
membrane_definition = .2; %Distance to membrane (in microns) to be considered a membrane spot

membrane_pixels = membrane_definition/.130; %.130 microns/pixel

percent_off_center = .1; % How Far off center axis (in percent of cell width) to be considered a middle pixel

rid = 1;
for ri = 1:length(cell_struct)
    if size(cell_struct(ri).Spots) == [0,0];
        
        continue
    end
    cell_region_struct(rid).Cell = ri;
    cell_region_struct(rid).Total_Spots = cell_struct(ri).Num_Spots;
    cell_region_struct(rid).Length = cell_struct(ri).Cell_Y_Axis;
    membrane = 0;
    middle = 0;
    cytoplasm = 0;
    
    cell_region_struct(rid).Whole_Volume = 0;
    for sk = 2:length(cell_struct(ri).Transformed_Boundaries{1,1})
        cell_region_struct(rid).Whole_Volume = cell_region_struct(rid).Whole_Volume + abs(cell_struct(ri).Transformed_Boundaries{1,1}(sk,1)-cell_struct(ri).Transformed_Boundaries{1,1}(sk-1,1))*(.5*(pi*cell_struct(ri).Transformed_Boundaries{1,1}(sk,2)^2));
    end
    cell_region_struct(rid).Whole_Volume = pi*(.5*cell_struct(ri).Cell_X_Axis)^2*(cell_struct(ri).Cell_Y_Axis-cell_struct(ri).Cell_X_Axis)+(4/3)*pi*(.5*cell_struct(ri).Cell_X_Axis)^3;
    d_pole = .5*cell_struct(ri).Cell_Y_Axis - (1-pole_definition)*.5*cell_struct(ri).Cell_Y_Axis;
    %pole_angle = acos((cell_struct(ri).Cell_X_Axis-d_pole)/cell_struct(ri).Cell_X_Axis);
    pole_volume = (1/3)*(pi*d_pole*(3*((cell_struct(ri).Cell_X_Axis)^2-(cell_struct(ri).Cell_X_Axis-d_pole)^2)+(d_pole)^2));
    cell_region_struct(rid).Pole_Volume = pole_volume;
    cell_region_struct(rid).Membrane_Volume = 2*((1-pole_definition)*cell_struct(ri).Cell_Y_Axis*(cell_struct(ri).Cell_X_Axis*membrane_pixels-membrane_pixels^2)) + pole_volume;
    %cell_region_struct(rid).Middle_Volume = middle_definition*cell_struct(ri).Cell_Y_Axis*pi*(.5*percent_off_center*cell_struct(ri).Cell_X_Axis)^2;
    cell_region_struct(rid).Middle_Volume = 0;
    for ji = 1:length(cell_struct(ri).Transformed_Dapi_Boundaries)
        for jid = 2:length(cell_struct(ri).Transformed_Dapi_Boundaries{ji,1})
            cell_region_struct(rid).Middle_Volume = cell_region_struct(rid).Middle_Volume + abs(cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid,1)-cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid-1,1))*(.5*(pi*cell_struct(ri).Transformed_Dapi_Boundaries{ji,1}(jid,2)^2));
        end
    end
    
    cell_region_struct(rid).Cytoplasm_Volume = max([0,cell_region_struct(rid).Whole_Volume - cell_region_struct(rid).Pole_Volume - cell_region_struct(rid).Membrane_Volume - cell_region_struct(rid).Middle_Volume]);
    
    for sri = 1:length(cell_struct(ri).Spots(:,1))
        from_pole = abs(cell_struct(ri).Spots(sri,3));
        if spot_struct(cell_struct(ri).Spots(sri,1)).Distance2Membrane <= membrane_pixels
            membrane = membrane + 1;
            cell_struct(ri).Spots(sri,5) = 1;
            continue
        elseif spot_struct(cell_struct(ri).Spots(sri,1)).Region == 4
            middle = middle + 1;
            cell_struct(ri).Spots(sri,5) = 4;
            continue
        elseif from_pole > (1-pole_definition)*.5*cell_struct(ri).Cell_Y_Axis
            membrane = membrane+1;
            cell_struct(ri).Spots(sri,5) = 1;
            continue
        else
            cytoplasm = cytoplasm + 1;
            cell_struct(ri).Spots(sri,5) = 2;
            continue
        end
    end
         
    cell_region_struct(rid).Pole_Spots = pole;
    cell_region_struct(rid).Membrane_Spots = membrane;
    cell_region_struct(rid).Middle_Spots = middle;
    cell_region_struct(rid).Cytoplasm_Spots = cytoplasm;
    rid = rid + 1;
    
end
    

