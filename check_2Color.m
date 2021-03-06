spot_struct = struct(field1,[],field2, [], field3, [], field4, [], field5, [], field6, [], field7, [],field8,[]); 

se = [1 1 1 ; 1 1 1 ; 1 1 1];
counted = [];

all_cells_good = 0;

for i = 1:max(new_mask(:))
    strcat(['Working on Cell ' num2str(i)])
    cell_struct(i).Cell = i;
    cell_struct(i).Center = [centers(i,2),centers(i,1)];
    cell_struct(i).Boundaries = bwboundaries(imdilate(new_mask==i,se));
    cell_struct(i).Dapi_Boundaries = bwboundaries(d4==i);
    cell_struct(i).Transformed_Dapi_Boundaries = cell_struct(i).Dapi_Boundaries;
    cell_struct(i).Cell_Angle = -ellipticity(i,2)*(pi/180);
    cell_struct(i).Cell_Y_Axis = ellipticity(i,3);
    cell_struct(i).Cell_X_Axis = ellipticity(i,4);
    cell_struct(i).Transformed_Boundaries = cell_struct(i).Boundaries;
    
    cell_struct(i).Transformed_Boundaries{1,1}(:,1) = cell_struct(i).Transformed_Boundaries{1,1}(:,1) - cell_struct(i).Center(2);
    cell_struct(i).Transformed_Boundaries{1,1}(:,2) = cell_struct(i).Transformed_Boundaries{1,1}(:,2) - cell_struct(i).Center(1);
    row_border = cell_struct(i).Transformed_Boundaries{1,1}(:,1);
    col_border = cell_struct(i).Transformed_Boundaries{1,1}(:,2);
    
    cell_struct(i).Transformed_Boundaries{1,1}(:,1) = (col_border*sin(cell_struct(i).Cell_Angle)+row_border*cos(cell_struct(i).Cell_Angle));
    cell_struct(i).Transformed_Boundaries{1,1}(:,2) = (col_border*cos(cell_struct(i).Cell_Angle)-row_border*sin(cell_struct(i).Cell_Angle));
    
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
    
    red_cell_spots = [];
    green_cell_spots = [];
    mask_i = imdilate(new_mask==i,se);
    dapi_mask_i = d4 == i;
    [cell_row, cell_col] = ind2sub(size(mask_i),find(mask_i));
    
    for g = 1:length(red_spot_coords)
        if sum(cell_row==red_spot_coords(g,2) & cell_col==red_spot_coords(g,3)) == 1
            red_cell_spots = [red_cell_spots red_spot_coords(g,1)];
            spot_struct(red_spot_coords(g,1)).Spot = [red_spot_coords(g,1), 1];
            spot_struct(red_spot_coords(g,1)).Coordinate = [red_x(red_spot_coords(g,1)),red_y(red_spot_coords(g,1)),red_z(red_spot_coords(g,1))];
            spot_struct(red_spot_coords(g,1)).Cell = i;
            spot_struct(red_spot_coords(g,1)).Region = 1;
        end
    end
    
    for gi = 1:length(green_spot_coords)-1
        if sum(cell_row==green_spot_coords(gi,2) & cell_col==green_spot_coords(gi,3)) == 1
            green_cell_spots = [green_cell_spots green_spot_coords(gi,1)];
            spot_struct(green_spot_coords(gi,1)+green_start).Spot = [green_spot_coords(gi,1), 3];
            spot_struct(green_spot_coords(gi,1)+green_start).Coordinate = [green_x(gi),green_y(gi),green_z(gi)];
            spot_struct(green_spot_coords(gi,1)+green_start).Region = 1;
        end
    end
    
    cell_struct(i).Num_Spots = [length(red_cell_spots),length(green_cell_spots)];
    red_cell_spot_coords = zeros(length(red_cell_spots),5);
    green_cell_spot_coords = zeros(length(green_cell_spots),5);
    
    qr = 1;
    for gs = red_cell_spots
        mask_i_temp = mask_i;
        
        for r = 1:3
            mask_i_temp = imerode(mask_i_temp,se);
            [cell_row_temp, cell_col_temp] = ind2sub(size(mask_i_temp),find(mask_i_temp));
            if sum(cell_row_temp==red_spot_coords_fixed(gs,2) & cell_col_temp==red_spot_coords_fixed(gs,3)) == 1
                spot_struct(gs).Region = spot_struct(gs).Region + 1;
            end
        end
        
        [cell_row_temp, cell_col_temp] = ind2sub(size(dapi_mask_i),find(dapi_mask_i));
        if sum(cell_row_temp==red_spot_coords_fixed(gs,2) & cell_col_temp==red_spot_coords_fixed(gs,3)) == 1
            spot_struct(gs).Region = 4;
        end
        
        
        spot_struct(gs).Transform_3D_Coordinate = spot_struct(gs).Coordinate;
        spot_struct(gs).Transform_3D_Coordinate(1) = spot_struct(gs).Transform_3D_Coordinate(1) - cell_struct(i).Center(2);
        spot_struct(gs).Transform_3D_Coordinate(2) = spot_struct(gs).Transform_3D_Coordinate(2) - cell_struct(i).Center(1);
        spot_row = spot_struct(gs).Transform_3D_Coordinate(1);
        spot_col = spot_struct(gs).Transform_3D_Coordinate(2);
        spot_struct(gs).Transform_3D_Coordinate(2) = (spot_col*sin(cell_struct(i).Cell_Angle)+spot_row*cos(cell_struct(i).Cell_Angle));
        spot_struct(gs).Transform_3D_Coordinate(1) = (spot_col*cos(cell_struct(i).Cell_Angle)-spot_row*sin(cell_struct(i).Cell_Angle));
        
        red_cell_spot_coords(qr,:) = [gs spot_struct(gs).Transform_3D_Coordinate(1) spot_struct(gs).Transform_3D_Coordinate(2) spot_struct(gs).Transform_3D_Coordinate(3) spot_struct(gs).Region];
        qr = qr+1;
    end
    
    qg = 1;
    for gr = green_cell_spots
        mask_i_temp = mask_i;
        
        for r = 1:3
            mask_i_temp = imerode(mask_i_temp,se);
            [cell_row_temp, cell_col_temp] = ind2sub(size(mask_i_temp),find(mask_i_temp));
            if sum(cell_row_temp==green_spot_coords_fixed(gr,2) & cell_col_temp==green_spot_coords_fixed(gr,3)) == 1
                spot_struct(gr+green_start).Region = spot_struct(gr+green_start).Region + 1;
            end
        end
        
        [cell_row_temp, cell_col_temp] = ind2sub(size(dapi_mask_i),find(dapi_mask_i));
        if sum(cell_row_temp==green_spot_coords_fixed(gr,2) & cell_col_temp==green_spot_coords_fixed(gr,3)) == 1
            spot_struct(gr+green_start).Region = 4;
        end
        
        
        spot_struct(gr+green_start).Transform_3D_Coordinate = spot_struct(gr+green_start).Coordinate;
        spot_struct(gr+green_start).Transform_3D_Coordinate(1) = spot_struct(gr+green_start).Transform_3D_Coordinate(1) - cell_struct(i).Center(2);
        spot_struct(gr+green_start).Transform_3D_Coordinate(2) = spot_struct(gr+green_start).Transform_3D_Coordinate(2) - cell_struct(i).Center(1);
        spot_row = spot_struct(gr+green_start).Transform_3D_Coordinate(1);
        spot_col = spot_struct(gr+green_start).Transform_3D_Coordinate(2);
        spot_struct(gr+green_start).Transform_3D_Coordinate(2) = (spot_col*sin(cell_struct(i).Cell_Angle)+spot_row*cos(cell_struct(i).Cell_Angle));
        spot_struct(gr+green_start).Transform_3D_Coordinate(1) = (spot_col*cos(cell_struct(i).Cell_Angle)-spot_row*sin(cell_struct(i).Cell_Angle));
        
        green_cell_spot_coords(qg,:) = [gr+green_start spot_struct(gr+green_start).Transform_3D_Coordinate(1) spot_struct(gr+green_start).Transform_3D_Coordinate(2) spot_struct(gr+green_start).Transform_3D_Coordinate(3) spot_struct(gr+green_start).Region];
        qg = qg+1;
    end
    
    red_cell_spot_coords_temp1 = red_cell_spot_coords;
    red_cell_spot_coords_temp2 = red_cell_spot_coords;
    green_cell_spot_coords_temp1 = green_cell_spot_coords;
    green_cell_spot_coords_temp2 = green_cell_spot_coords;
    
    big_axis = max(red_cell_spot_coords(:,3))-min(red_cell_spot_coords(:,3));
    small_axis = min(red_cell_spot_coords(:,2))-min(red_cell_spot_coords(:,2));
    
    if big_axis < small_axis
        red_cell_spot_coords_temp1(:,2) = red_cell_spot_coords_temp2(:,3);
        red_cell_spot_coords_temp1(:,3) = red_cell_spot_coords_temp2(:,2);
        green_cell_spot_coords_temp1(:,2) = green_cell_spot_coords_temp2(:,3);
        green_cell_spot_coords_temp1(:,3) = green_cell_spot_coords_temp2(:,2);
    end
    
    red_cell_spot_coords = red_cell_spot_coords_temp1;
    green_cell_spot_coords = green_cell_spot_coords_temp1;
    
    scatter(red_cell_spot_coords(:,2),red_cell_spot_coords(:,4),100,'filled','r');
    hold on
    scatter(green_cell_spot_coords(:,2),green_cell_spot_coords(:,4),100,'filled','g');
    grid on;title(strcat(['Cell ' num2str(i)]));xlabel('X Axis');ylabel('Z Axis')
    continuetranslation=1;
    discard_cell = 0;
    while continuetranslation==1 && all_cells_good == 0
        prompt={'X_Translation(- = Left, + = Right):','Y_Translation(+ = Up, - = Down):','Continue the translation(1 == Yes && 0== NO) ?','Discard Cell (1 == Yes) ?','Skip Rest of Cells(1 == Yes) ?'};   % A box will take in the values for the X/Ytranslation
        cell_title=strcat(['Cell ' num2str(i) ' Check']);                             % The title of the box
        answer=inputdlg(prompt,cell_title);
        Xtranslation = str2num(answer{1});
        Ztranslation = str2num(answer{2});
        continuetranslation = str2num(answer{3});
        discard_cell = str2num(answer{4});
        if discard_cell == 1
            break
        end
        all_cells_good = str2num(answer{5});
        if isempty(answer{5})
            all_cells_good = 0;
        end
        red_cell_spot_coords(:,2)=red_cell_spot_coords(:,2)+Xtranslation;                                % Translated
        red_cell_spot_coords(:,4)=red_cell_spot_coords(:,4)+Ztranslation;
        green_cell_spot_coords(:,2)=green_cell_spot_coords(:,2)+Xtranslation;                                % Translated
        green_cell_spot_coords(:,4)=green_cell_spot_coords(:,4)+Ztranslation;
        scatter(red_cell_spot_coords(:,2),red_cell_spot_coords(:,4),100,'filled','r');
        hold on
        scatter(green_cell_spot_coords(:,2),green_cell_spot_coords(:,4),100,'filled','g');
        grid on;title(strcat(['Cell ' num2str(i)]));xlabel('X Axis');ylabel('Z Axis')
    end
    
    if discard_cell == 1
        for gs = red_cell_spots
            
            spot_struct(gs).Cell = [];
            spot_struct(gs).Region = [];
            spot_struct(gs).Transform_3D_Coordinate = [];
            spot_struct(gs).Collapsed_2D_Coordinate = [];
            red_spot_coords(red_spot_coords(:,1)==gs,:) = [];
            
            
        end
        
        for gr = green_cell_spots
            
            spot_struct(gr+green_start).Cell = [];
            spot_struct(gr+green_start).Region = [];
            spot_struct(gr+green_start).Transform_3D_Coordinate = [];
            spot_struct(gr+green_start).Collapsed_2D_Coordinate = [];
            green_spot_coords(green_spot_coords(:,1)==gr,:) = [];
            
        end
        
            continue
    end
    
    cell_struct(i).Red_Spots = red_cell_spot_coords;
    cell_struct(i).Green_Spots = green_cell_spot_coords;
    
    for qt = 1:length(red_cell_spot_coords)
        
        spot_struct(red_cell_spot_coords(qt,1)).Transform_3D_Coordinate(2) = red_cell_spot_coords(qt,3);
        spot_struct(red_cell_spot_coords(qt,1)).Transform_3D_Coordinate(1) = red_cell_spot_coords(qt,2);
        spot_struct(red_cell_spot_coords(qt,1)).Transform_3D_Coordinate(3) = red_cell_spot_coords(qt,4);
        spot_struct(red_cell_spot_coords(qt,1)).Collapsed_2D_Coordinate = spot_struct(red_cell_spot_coords(qt,1)).Transform_3D_Coordinate;
        spot_struct(red_cell_spot_coords(qt,1)).Collapsed_2D_Coordinate(1) = sqrt((double(spot_struct(red_cell_spot_coords(qt,1)).Transform_3D_Coordinate(1)))^2+(double(spot_struct(red_cell_spot_coords(qt,1)).Transform_3D_Coordinate(3)))^2);
        spot_struct(red_cell_spot_coords(qt,1)).Collapsed_2D_Coordinate(3) = 0;
        spot_struct(red_cell_spot_coords(qt,1)).Distance2Center = spot_struct(red_cell_spot_coords(qt,1)).Coordinate;
        spot_struct(red_cell_spot_coords(qt,1)).Distance2Center(1) = abs(.5*spot_struct(red_cell_spot_coords(qt,1)).Transform_3D_Coordinate(1)/cell_struct(i).Cell_X_Axis);
        spot_struct(red_cell_spot_coords(qt,1)).Distance2Center(2) = abs(.5*spot_struct(red_cell_spot_coords(qt,1)).Transform_3D_Coordinate(2)/cell_struct(i).Cell_Y_Axis);
        spot_struct(red_cell_spot_coords(qt,1)).Distance2Center(3) = abs(.5*spot_struct(red_cell_spot_coords(qt,1)).Transform_3D_Coordinate(3)/cell_struct(i).Cell_X_Axis);
        spot_struct(red_cell_spot_coords(qt,1)).Distance2Membrane = Distance2Edge(cell_struct(i).Transformed_Boundaries,[spot_struct(red_cell_spot_coords(qt,1)).Collapsed_2D_Coordinate(2),spot_struct(red_cell_spot_coords(qt,1)).Collapsed_2D_Coordinate(1)]);
        red_spot_coords(red_spot_coords(:,1)==red_cell_spot_coords(qt,1),:) = [];
        
    end
    
    for gt = 1:length(green_cell_spot_coords)
        
        spot_struct(green_cell_spot_coords(gt,1)).Transform_3D_Coordinate(2) = green_cell_spot_coords(gt,3);
        spot_struct(green_cell_spot_coords(gt,1)).Transform_3D_Coordinate(1) = green_cell_spot_coords(gt,2);
        spot_struct(green_cell_spot_coords(gt,1)).Transform_3D_Coordinate(3) = green_cell_spot_coords(gt,4);
        spot_struct(green_cell_spot_coords(gt,1)).Collapsed_2D_Coordinate = spot_struct(green_cell_spot_coords(gt,1)).Transform_3D_Coordinate;
        spot_struct(green_cell_spot_coords(gt,1)).Collapsed_2D_Coordinate(1) = sqrt((double(spot_struct(green_cell_spot_coords(gt,1)).Transform_3D_Coordinate(1)))^2+(double(spot_struct(green_cell_spot_coords(gt,1)).Transform_3D_Coordinate(3)))^2);
        spot_struct(green_cell_spot_coords(gt,1)).Collapsed_2D_Coordinate(3) = 0;
        spot_struct(green_cell_spot_coords(gt,1)).Distance2Center = spot_struct(green_cell_spot_coords(gt,1)).Coordinate;
        spot_struct(green_cell_spot_coords(gt,1)).Distance2Center(1) = abs(.5*spot_struct(green_cell_spot_coords(gt,1)).Transform_3D_Coordinate(1)/cell_struct(i).Cell_X_Axis);
        spot_struct(green_cell_spot_coords(gt,1)).Distance2Center(2) = abs(.5*spot_struct(green_cell_spot_coords(gt,1)).Transform_3D_Coordinate(2)/cell_struct(i).Cell_Y_Axis);
        spot_struct(green_cell_spot_coords(gt,1)).Distance2Center(3) = abs(.5*spot_struct(green_cell_spot_coords(gt,1)).Transform_3D_Coordinate(3)/cell_struct(i).Cell_X_Axis);
        spot_struct(green_cell_spot_coords(gt,1)).Distance2Membrane = Distance2Edge(cell_struct(i).Transformed_Boundaries,[spot_struct(green_cell_spot_coords(gt,1)).Collapsed_2D_Coordinate(2),spot_struct(green_cell_spot_coords(gt,1)).Collapsed_2D_Coordinate(1)]);
        green_spot_coords(green_spot_coords(:,1)==green_cell_spot_coords(gt,1),:) = [];
        
    end
    
    
end
    
 