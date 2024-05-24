function [bend] = findbend(Panel, Node)
bend = zeros(length(Panel),4);
for i = 1:length(Panel)
    if numel(Panel{i}) == 4
        L1 = norm((Node(Panel{i}(1),:)-Node(Panel{i}(3),:)));
        L2 = norm((Node(Panel{i}(4),:)-Node(Panel{i}(2),:)));
        if L1>L2, lclbend = [2,4,1,3];
        else lclbend = [1,3,2,4]; end;
        bend(i,:) = Panel{i}(lclbend);
    end
end;
bend(sum(bend,2)==0,:)=[];  