X = 256;
Y = 256;

% for i = length(spot_struct):-1:1
%     
%     
%     if isempty(spot_struct(i).Distance2Center)
%         spot_struct(i) = [];
%     end
%     
% end

x1 = [];
y1 = [];
x2 = [];
y2 = [];
x3 = [];
y3 = [];
x4 = [];
y4 = [];


for j = 1:length(spot_struct) 
    if spot_struct(j).Region == 1
        x1(j) = int32(spot_struct(j).Coordinate(1));
        y1(j) = int32(spot_struct(j).Coordinate(2));
    elseif spot_struct(j).Region == 2
        x2(j) = int32(spot_struct(j).Coordinate(1));
        y2(j) = int32(spot_struct(j).Coordinate(2));
    elseif spot_struct(j).Region == 3
        x3(j) = int32(spot_struct(j).Coordinate(1));
        y3(j) = int32(spot_struct(j).Coordinate(2));
    elseif spot_struct(j).Region == 4
        x4(j) = int32(spot_struct(j).Coordinate(1));
        y4(j) = int32(spot_struct(j).Coordinate(2));
    end
end

x1(x1<1) = 1;
y1(y1<1) = 1;
x1(x1>=X) = X;
y1(y1>=Y) = Y;

x2(x2<1) = 1;
y2(y2<1) = 1;
x2(x2>=X) = X;
y2(y2>=Y) = Y;

x3(x3<1) = 1;
y3(y3<1) = 1;
x3(x3>=X) = X;
y3(y3>=Y) = Y;

x4(x4<1) = 1;
y4(y4<1) = 1;
x4(x4>=X) = X;
y4(y4>=Y) = Y;


figure(1)
scatter(y1,x1,10,'r','filled');
figure(2)
scatter(y2,x2,10,'r','filled');
figure(3)
scatter(y3,x3,10,'r','filled');
figure(4)
scatter(y4,x4,10,'r','filled');


