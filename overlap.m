function  overlap( imagefilename,x,y,Color)
image = imread(imagefilename);                                            % Put the reference image based on which you want to select the points
imshow(image);
set(gcf, 'Position', get(0,'Screensize'));                                % Maximizing the image screenwindow
hold on     
DifferentColorCodes=unique(Color);
cc=hsv(numel(DifferentColorCodes));
for m=1:numel(DifferentColorCodes) 
    valuesIndexes=find(Color==DifferentColorCodes(m));                    % Finding Indices which have the particular ColorCode
    x1=x(valuesIndexes);                                                  % Choosing the x,y values which belong to the particular indices (i.e to a particular color code)
    y1=y(valuesIndexes);
plot(x1,y1,'.','color',cc(m,:),'LineWidth',0.005);                        % Plotting the points and cc(m:,) is just giving the points different colors for different set of points ie ( belonging to a paricular colorcode)
end
end


