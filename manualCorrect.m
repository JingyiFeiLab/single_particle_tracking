function newObjects = manualCorrect(fig,objects);

fig;
p = get(fig, 'Position');

%non_mask = 2.*non_mask;
while 0<1
    [y,x,button] = ginput(1);
    if isempty(button) % Enter
        break;
%     elseif button == 112 % p for false positive selevtion
%         %[false_y, false_x] = ginput();
%         false_y = y;
%         false_x = x;
%         false_y = int64(false_y);
%         false_x = int64(false_x);
%         false_positives = zeros(length(false_x),1);
%         
%         for i = 1:length(false_x)
%             false_positives(i) = objects(false_x(i),false_y(i));
%             mask(objects == false_positives(i)) = 0 ;
%             non_mask(objects == false_positives(i)) = 1;
%         end
%         
%         IR = original;
%         IG = original;
%         IB = original;
%         
%         IR(bwperim(all)) = 0;
%         IG(bwperim(all)) = 0;
%         IB(bwperim(all)) = 0;
%         
%         IR(bwperim(non_mask)) = 100;
%         IB(bwperim(mask)) = 100;
%         
%         IRGB = cat(3,IR,IG,IB);
%         imshow(IRGB);
%         set(figure,'Position',p)
%         set(gca,'Units','normalized', 'Position' ,[0.05 0.05 .95 .85])
%         title({'Press "p" to switch positive selections to negative , Press "n" to switch negative to positive' ...
%                                           'Press "d" to delete selections , Press "s" to manually split cells'...
%                                           'Press "w"  to manually draw border (USE CROSSHAIRS TO PICK CELL)'})
%                                       
%         str = {'Positive Selection', 'Negative Selection'};
%         t = text([2 2],[4 8],str);
%         t(1).Color = 'blue'; t(1).FontSize = 14 ;
%         t(2).Color = 'red'; t(2).FontSize = 14 ;                              
%                                       
%         if exist('new_borders') == 0
%             varargout{1} = [];
%         end
        
%     elseif button == 110; % n for false negatives
%         
%         %[false_neg_y false_neg_x] = ginput();
%         false_neg_y = y;
%         false_neg_x = x;
%         false_neg_y = int64(false_neg_y);
%         false_neg_x = int64(false_neg_x);
%         false_negatives = zeros(length(false_neg_x),1);
%         
%         for i = 1:length(false_neg_x)
%             false_negatives(i) = objects(false_neg_x(i),false_neg_y(i));
%             mask = mask + (objects == false_negatives(i)) ;
%             non_mask = non_mask - (objects == false_negatives(i));
%         end
%         
%         IR = original;
%         IG = original;
%         IB = original;
%         
%         IR(bwperim(all)) = 0;
%         IG(bwperim(all)) = 0;
%         IB(bwperim(all)) = 0;
%         
%         IB(bwperim(mask)) = 100 ;
%         IR(bwperim(non_mask)) = 100 ;
%         
%         
%         IRGB = cat(3,IR,IG,IB);
%         imshow(IRGB);title({'Press "p" to switch positive selections to negative , Press "n" to switch negative to positive' ...
%                                           'Press "d" to delete selections , Press "s" to manually split cells'...
%                                           'Press "w"  to manually draw border (USE CROSSHAIRS TO PICK CELL)'})
%                                       
%         set(figure,'Position',p)
%         set(gca,'Units','normalized', 'Position' ,[0.05 0.05 .95 .85])
%         
%         str = {'Positive Selection', 'Negative Selection'};
%         t = text([2 2],[4 8],str);
%         t(1).Color = 'blue'; t(1).FontSize = 14 ;
%         t(2).Color = 'red'; t(2).FontSize = 14 ;         
%         
%         if exist('new_borders') == 0
%             varargout{1} = [];
%         end
    
    elseif button == 100 % d for delete
        
        %[delete_y, delete_x] = ginput();
        delete_y = y;
        delete_x = x;
        delete_y = int64(delete_y);
        delete_x = int64(delete_x);
        deletions = zeros(length(delete_y),1);
        
        for i = 1:length(delete_y)
            deletions(i) = objects(delete_x(i),delete_y(i));
            objects(objects == deletions(i)) = 0;
        end
        
        IR = objects;
        IG = objects;
        IB = objects;
        
        IR(bwperim(objects)) = 100;
        IG(bwperim(objects)) = 100;
        IB(bwperim(objects)) = 100;
        
        IRGB = cat(3,IR,IG,IB);
        imshow(IRGB);title({'Press "d" to delete selections , Press "s" to manually split cells'...
                            'Press "w"  to manually draw border (USE CROSSHAIRS TO PICK CELL)' })
                         
        set(fig,'Position',p)
        set(gca,'Units','normalized', 'Position' ,[0.05 0.05 .95 .85])
       
        if exist('new_borders') == 0
            varargout{1} = [];
        end
        

            
    elseif button == 119 % w for manual draw border
        new_borders = newBorder(fig);
        objects = newObject(objects,new_borders);
        
        
        IR = objects;
        IG = objects;
        IB = objects;
        
        IR(bwperim(objects)) = 100;
        IG(bwperim(objects)) = 100;
        IB(bwperim(objects)) = 100;
        
        
        IRGB = cat(3,IR,IG,IB);
        imshow(IRGB);title({'Press "d" to delete selections , Press "s" to manually split cells'...
                            'Press "w"  to manually draw border (USE CROSSHAIRS TO PICK CELL)' })
                                      
        set(fig,'Position',p)
        set(gca,'Units','normalized', 'Position' ,[0.05 0.05 .95 .85])
        
        
    elseif button == 115 % s for split
        
        %objects = bwlabel(objects,4);
        objects1 = objects;
        h = imfreehand('Closed',false);
        x = int64(getPosition(h));
        new_line.line = x;
        
        for i = 1:length(new_line(1).line)
            for j = 1:length(new_line(1).line)
                objects(new_line(1).line(j,2),new_line(1).line(j,1)) = 0;
            end
%             new_object = imfill(new_object,'holes');
        end
        
        objects = bwlabel(objects);

%         split_start = x(1,:);
%         split_end = x(length(x),:);
%         bad_object = max(max(objects(x(:,2),x(:,1))));
%         bound = bwboundaries(objects == bad_object);
%         objects = objects - (objects == bad_object);
%         bad_object = (objects1 == bad_object);
%         
%         
%         edge_pic1 = zeros(length(bound{1,1}(:,1)),1);
%         edge_pic2 = zeros(length(bound{1,1}(:,1)),1);
%         
%         for jd = 1:length(bound{1,1}(:,1))
%             edge_pic1(jd) = distance(bound{1,1}(jd,:),[split_start(2),split_start(1)]);
%         end
%         
%         [~,e1] = min(edge_pic1);
%         
%         for kd = 1:length(bound{1,1}(:,1))
%             edge_pic2(kd) = distance(bound{1,1}(kd,:),[split_end(2),split_end(1)]);
%         end
%         
%         [~,e2] = min(edge_pic2);
%         
%         xv1 = [bound{1,1}(e1,2),split_start(2)];
%         xv2 = [bound{1,1}(e2,2),split_end(2)];
%         yv1 = [bound{1,1}(e1,1),split_start(1)];
%         yv2 = [bound{1,1}(e2,1),split_start(2)];
%         
%         line1 = makeLine(xv1,yv1);
%         line2 = makeLine(xv2,yv2);
%         
%         bad_object(sub2ind(size(objects),line1(:,2),line1(:,1))) = 0;
%         bad_object(sub2ind(size(objects),line2(:,2),line2(:,1))) = 0;
%         bad_object(sub2ind(size(objects),x(:,2),x(:,1))) = 0;
%         bad_object = bwlabel(bad_object,4);
%         
%        
%         max_id = max(max(objects));
%         objects = objects + ((max_id+bad_object).*im2bw(bad_object,1));
        IR = objects;
        IG = objects;
        IB = objects;
        
    
        IR(bwperim(objects)) = 100;
        IG(bwperim(objects)) = 100;
        IB(bwperim(objects)) = 100;
        
        IRGB = cat(3,IR,IG,IB);
        imshow(IRGB);title({'Press "d" to delete selections , Press "s" to manually split cells'...
                            'Press "w"  to manually draw border (USE CROSSHAIRS TO PICK CELL)'})
        
        set(fig,'Position',p)
        set(gca,'Units','normalized', 'Position' ,[0.05 0.05 .95 .85])
        
        
        
        
        if exist('new_borders') == 0
            varargout{1} = [];
        end
        
    end
    
end



indices = unique(objects(:));
index = 0;

for r = 2:length(indices)
    index = index + 1 ;
    objects(objects == indices(r)) = index;
end

newObjects = objects;



end