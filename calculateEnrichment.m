membrane_enrichment = [];
cytoplasm_enrichment = [];
pole_enrichment = [];
middle_enrichment = [];
non_middle_enrichment = [];

membrane_density = [];
cytoplasm_density = [];
pole_density = [];
middle_density = [];
non_middle_density = [];

cell_num = [];
outliers = [];

for i = 1:length(cell_region_struct)
    
    if ismember(cell_region_struct(i).Cell,outliers);
        continue
    end
    
    if cell_region_struct(i).Middle_Volume == 0
        continue
    else
        middle_enrichment = [middle_enrichment (cell_region_struct(i).Middle_Spots/cell_region_struct(i).Total_Spots)/(cell_region_struct(i).Middle_Volume/cell_region_struct(i).Whole_Volume)];
        middle_density = [middle_density (cell_region_struct(i).Middle_Spots/cell_region_struct(i).Middle_Volume)];
    end
    
    if cell_region_struct(i).Membrane_Volume == 0
        membrane_enrichment = [membrane_enrichment 0];
        membrane_density = [membrane_density 0];
    else
        membrane_enrichment = [membrane_enrichment (cell_region_struct(i).Membrane_Spots/cell_region_struct(i).Total_Spots)/(cell_region_struct(i).Membrane_Volume/cell_region_struct(i).Whole_Volume)];
        membrane_density = [membrane_density (cell_region_struct(i).Membrane_Spots/cell_region_struct(i).Membrane_Volume)];
    end
    
    if cell_region_struct(i).Cytoplasm_Volume == 0
        cytoplasm_enrichment = [cytoplasm_enrichment 0];
        cytoplasm_density = [cytoplasm_density 0];
    else
        cytoplasm_enrichment = [cytoplasm_enrichment (cell_region_struct(i).Cytoplasm_Spots/cell_region_struct(i).Total_Spots)/(cell_region_struct(i).Cytoplasm_Volume/cell_region_struct(i).Whole_Volume)];
        cytoplasm_density = [cytoplasm_density (cell_region_struct(i).Cytoplasm_Spots/cell_region_struct(i).Cytoplasm_Volume)];
    end
    
    if cell_region_struct(i).Pole_Volume == 0
        pole_enrichment = [pole_enrichment 0];
        pole_density = [pole_density 0];
    else
        pole_enrichment = [pole_enrichment (cell_region_struct(i).Pole_Spots/cell_region_struct(i).Total_Spots)/(cell_region_struct(i).Pole_Volume/cell_region_struct(i).Whole_Volume)];
        pole_density = [pole_density (cell_region_struct(i).Pole_Spots/cell_region_struct(i).Pole_Volume)];
    end
    
    non_middle_enrichment = [non_middle_enrichment ((cell_region_struct(i).Membrane_Spots+cell_region_struct(i).Cytoplasm_Spots+cell_region_struct(i).Pole_Spots)/cell_region_struct(i).Total_Spots)/((cell_region_struct(i).Membrane_Volume+cell_region_struct(i).Cytoplasm_Volume+cell_region_struct(i).Pole_Volume)/cell_region_struct(i).Whole_Volume)];
    non_middle_density = [non_middle_density ((cell_region_struct(i).Membrane_Spots+cell_region_struct(i).Cytoplasm_Spots+cell_region_struct(i).Pole_Spots)/cell_region_struct(i).Membrane_Volume)];
    
    cell_num = [cell_num cell_region_struct(i).Cell];
end

cell_num(isnan(membrane_enrichment)) = [];

membrane_enrichment(isnan(membrane_enrichment)) = [];

cytoplasm_enrichment(isnan(cytoplasm_enrichment)) = [];

pole_enrichment(isnan(pole_enrichment)) = [];

middle_enrichment(isnan(middle_enrichment)) = [];

non_middle_enrichment(isnan(non_middle_enrichment)) = [];

membrane_density(isnan(membrane_density)) = [];

cytoplasm_density(isnan(cytoplasm_density)) = [];

pole_density(isnan(pole_density)) = [];

middle_density(isnan(middle_density)) = [];

non_middle_density(isnan(non_middle_density)) = [];

enrichment = [membrane_enrichment;cytoplasm_enrichment;pole_enrichment;middle_enrichment;cell_num];
mean_enrichment = [mean(membrane_enrichment) std(membrane_enrichment);mean(cytoplasm_enrichment) std(cytoplasm_enrichment);mean(pole_enrichment) std(pole_enrichment);mean(middle_enrichment) std(middle_enrichment)];

just_middle_enrichment = [non_middle_enrichment;middle_enrichment;cell_num];
just_middle_mean_enrichment = [mean(non_middle_enrichment) std(non_middle_enrichment);mean(middle_enrichment) std(middle_enrichment)];


save(filename)

% density = [membrane_density;cytoplasm_density;pole_density;middle_density;cell_num];
% mean_density = [mean(membrane_density) std(membrane_density);mean(cytoplasm_density) std(cytoplasm_density);mean(pole_density) std(pole_density);mean(middle_density) std(middle_density)];

