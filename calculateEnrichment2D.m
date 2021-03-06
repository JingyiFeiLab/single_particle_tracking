membrane_enrichment = [];
cytoplasm_enrichment = [];
pole_enrichment = [];
middle_enrichment = [];
not_middle_enrichment = [];

enrichment_2D = [];
mean_enrichment_2D = [];
just_middle_enrichment_2D = [];
just_middle_mean_enrichment_2D = [];

membrane_density = [];
cytoplasm_density = [];
pole_density = [];
middle_density = [];
not_middle_density = [];

cell_num = [];
outliers = [6,13,21,26,45,48,43];

for i = 1:length(cell_region_struct_2d)
    
    if ismember(cell_region_struct_2d(i).Cell,outliers);
        continue
    end
    
    if cell_region_struct_2d(i).Middle_Area == 0
        continue
    else
        middle_enrichment = [middle_enrichment (cell_region_struct_2d(i).Middle_Spots/cell_region_struct_2d(i).Total_Spots)/(cell_region_struct_2d(i).Middle_Area/cell_region_struct_2d(i).Whole_Area)];
        middle_density = [middle_density (cell_region_struct_2d(i).Middle_Spots/cell_region_struct_2d(i).Middle_Area)];
    end
    
    if cell_region_struct_2d(i).Membrane_Area == 0
        membrane_enrichment = [membrane_enrichment 0];
        membrane_density = [membrane_density 0];
    else
        membrane_enrichment = [membrane_enrichment (cell_region_struct_2d(i).Membrane_Spots/cell_region_struct_2d(i).Total_Spots)/(cell_region_struct_2d(i).Membrane_Area/cell_region_struct_2d(i).Whole_Area)];
        membrane_density = [membrane_density (cell_region_struct_2d(i).Membrane_Spots/cell_region_struct_2d(i).Membrane_Area)];
    end
    
    if cell_region_struct_2d(i).Cytoplasm_Area == 0
        cytoplasm_enrichment = [cytoplasm_enrichment 0];
        cytoplasm_density = [cytoplasm_density 0];
    else
        cytoplasm_enrichment = [cytoplasm_enrichment (cell_region_struct_2d(i).Cytoplasm_Spots/cell_region_struct_2d(i).Total_Spots)/(cell_region_struct_2d(i).Cytoplasm_Area/cell_region_struct_2d(i).Whole_Area)];
        cytoplasm_density = [cytoplasm_density (cell_region_struct_2d(i).Cytoplasm_Spots/cell_region_struct_2d(i).Cytoplasm_Area)];
    end
    
    if cell_region_struct_2d(i).Pole_Area == 0
        pole_enrichment = [pole_enrichment 0];
        pole_density = [pole_density 0];
    else
        pole_enrichment = [pole_enrichment (cell_region_struct_2d(i).Pole_Spots/cell_region_struct_2d(i).Total_Spots)/(cell_region_struct_2d(i).Pole_Area/cell_region_struct_2d(i).Whole_Area)];
        pole_density = [pole_density (cell_region_struct_2d(i).Pole_Spots/cell_region_struct_2d(i).Pole_Area)];
    end

    not_middle_enrichment = [not_middle_enrichment ((cell_region_struct_2d(i).Pole_Spots+cell_region_struct_2d(i).Cytoplasm_Spots+cell_region_struct_2d(i).Membrane_Spots)/cell_region_struct_2d(i).Total_Spots)/((cell_region_struct_2d(i).Pole_Area+cell_region_struct_2d(i).Cytoplasm_Area+cell_region_struct_2d(i).Membrane_Area)/cell_region_struct_2d(i).Whole_Area)];
    not_middle_density = [not_middle_density ((cell_region_struct_2d(i).Pole_Spots+cell_region_struct_2d(i).Cytoplasm_Spots+cell_region_struct_2d(i).Membrane_Spots)/(cell_region_struct_2d(i).Pole_Area+cell_region_struct_2d(i).Cytoplasm_Area+cell_region_struct_2d(i).Membrane_Area))];
    
    cell_num = [cell_num cell_region_struct_2d(i).Cell];
end

cell_num(isnan(membrane_enrichment)) = [];

membrane_enrichment(isnan(membrane_enrichment)) = [];

cytoplasm_enrichment(isnan(cytoplasm_enrichment)) = [];

pole_enrichment(isnan(pole_enrichment)) = [];

middle_enrichment(isnan(middle_enrichment)) = [];

not_middle_enrichment(isnan(not_middle_enrichment)) = [];

membrane_density(isnan(membrane_density)) = [];

cytoplasm_density(isnan(cytoplasm_density)) = [];

pole_density(isnan(pole_density)) = [];

middle_density(isnan(middle_density)) = [];

not_middle_density(isnan(not_middle_density)) = [];

enrichment_2D = [membrane_enrichment;cytoplasm_enrichment;pole_enrichment;middle_enrichment;cell_num];
mean_enrichment_2D = [mean(membrane_enrichment) std(membrane_enrichment);mean(cytoplasm_enrichment) std(cytoplasm_enrichment);mean(pole_enrichment) std(pole_enrichment);mean(middle_enrichment) std(middle_enrichment)];

just_middle_enrichment_2D = [not_middle_enrichment;middle_enrichment;cell_num];
just_middle_mean_enrichment_2D = [mean(not_middle_enrichment) std(not_middle_enrichment);mean(middle_enrichment) std(middle_enrichment)];

density2D = [membrane_density;cytoplasm_density;pole_density;middle_density;cell_num];
mean_density2D = [mean(membrane_density) std(membrane_density);mean(cytoplasm_density) std(cytoplasm_density);mean(pole_density) std(pole_density);mean(middle_density) std(middle_density)];

save(filename)