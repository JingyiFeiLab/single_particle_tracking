%close all
one_cell_id = 9;

figure(1)
plot(cell_struct(one_cell_id).Transformed_Boundaries{1,1}(:,2),cell_struct(one_cell_id).Transformed_Boundaries{1,1}(:,1),'k','LineWidth',3)
hold on
% for i = 1:length(cell_struct(one_cell_id).Transformed_Dapi_Boundaries)
%     plot(cell_struct(one_cell_id).Transformed_Dapi_Boundaries{i,1}(:,2),cell_struct(one_cell_id).Transformed_Dapi_Boundaries{i,1}(:,1),'m','LineWidth',3)
%     hold on
% end

membrane_boundaries = zeros(length(cell_struct(one_cell_id).Transformed_Boundaries{1,1}),2);


for sk = 1:length(cell_struct(one_cell_id).Transformed_Boundaries{1,1})
    point_angle = atan(cell_struct(one_cell_id).Transformed_Boundaries{1,1}(sk,1)/cell_struct(one_cell_id).Transformed_Boundaries{1,1}(sk,2));
    h = sqrt((cell_struct(one_cell_id).Transformed_Boundaries{1,1}(sk,1))^2+(cell_struct(one_cell_id).Transformed_Boundaries{1,1}(sk,2))^2);
    h_prime = h - membrane_pixels;
    sign_x = sign(cell_struct(one_cell_id).Transformed_Boundaries{1,1}(sk,1));
    sign_y = sign(cell_struct(one_cell_id).Transformed_Boundaries{1,1}(sk,2));
    x_prime = abs(h_prime*sin(point_angle));
    y_prime = abs(h_prime*cos(point_angle));
    membrane_boundaries(sk,1) = sign_y*y_prime;
    membrane_boundaries(sk,2) = sign_x*x_prime;
end

plot(membrane_boundaries(:,1),membrane_boundaries(:,2),'g','LineWidth',3)
hold on

pole_x = -1*(cell_struct(one_cell_id).Cell_X_Axis-1):(cell_struct(one_cell_id).Cell_X_Axis-1);
d_pole = .5*cell_struct(one_cell_id).Cell_Y_Axis - (1-pole_definition)*.5*cell_struct(one_cell_id).Cell_Y_Axis;
pole_y = (max(cell_struct(one_cell_id).Transformed_Boundaries{1,1}(:,1))-d_pole)*ones(1,length(-1*(cell_struct(one_cell_id).Cell_X_Axis-1):(cell_struct(one_cell_id).Cell_X_Axis-1)));
neg_pole_y = (min(cell_struct(one_cell_id).Transformed_Boundaries{1,1}(:,1))+d_pole)*ones(1,length(-1*(cell_struct(one_cell_id).Cell_X_Axis-1):(cell_struct(one_cell_id).Cell_X_Axis-1)));

plot(pole_x,pole_y,'b','LineWidth',2)
hold on
plot(pole_x,neg_pole_y,'b','LineWidth',2)

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 cell_struct(one_cell_id).Cell_X_Axis cell_struct(one_cell_id).Cell_Y_Axis];

if sum(cell_region==3) ~= 0
        figure(3);gscatter(cell_center_x,cell_center_z,cell_region,'kgbm','.',10,'filled');grid on
    else
        figure(3);gscatter(cell_center_x,cell_center_z,cell_region,'kgm','.',10,'filled');grid on
    end
    title(strcat(['Cell ' num2str(one_cell_id) ' XZ Spot Projection']),'FontSize',24)
    ylabel('Vertical Axis (Z)','FontSize',24)
    xlabel('Short Axis (X)','FontSize',24)