function [B, L] = dirc3d(Node,Ele)
Ne = size(Ele,1); Nn = size(Node,1);
D = [Node(Ele(:,2),1)-Node(Ele(:,1),1), Node(Ele(:,2),2)-Node(Ele(:,1),2),...
    Node(Ele(:,2),3)-Node(Ele(:,1),3)];
L = sqrt(D(:,1).^2+D(:,2).^2+D(:,3).^2);
D = [D(:,1)./L D(:,2)./L D(:,3)./L];
B = sparse(repmat((1:Ne)',1,6),[3*Ele(:,1)-2 3*Ele(:,1)-1 3*Ele(:,1),...
           3*Ele(:,2)-2 3*Ele(:,2)-1 3*Ele(:,2)],[D -D],Ne,3*Nn);
B = -B;
