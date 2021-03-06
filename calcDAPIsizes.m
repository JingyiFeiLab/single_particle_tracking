dapi_volumes_scaled = [];
dapi_sizes_2d = [];
for i = 1:length(cell_region_struct)
    if cell_region_struct(i).Middle_Volume == 0 || cell_region_struct(i).Whole_Volume == 0
        continue
    end
    dapi_volumes_scaled = [dapi_volumes_scaled cell_region_struct(i).Middle_Volume/cell_region_struct(i).Whole_Volume];
    dapi_sizes_2d = [dapi_sizes_2d cellArea2(d4,cell_region_struct(i).Cell)];

end
