function [ p0,index ] = T6_pretreatment( p,t )
edge = [t(:,[1 2]);t(:,[2 3]);t(:,[1 3])];
edge = unique(sort(edge,2),'rows');id = 1:size(edge,1);
facecen = [(p(edge(:,1),1) + p(edge(:,2),1))/2 (p(edge(:,1),2) + p(edge(:,2),2))/2];
midnode = size(p,1)+1:size(p,1)+size(edge,1);
for e = 1:size(t,1)
    edgee = [unique(t(e,[1 2]));unique(t(e,[2 3]));unique(t(e,[1 3]))];
    id1 = find(sum([edgee(1,1) == edge(:,1), edgee(1,2) == edge(:,2)],2) == 2);
    id2 = find(sum([edgee(2,1) == edge(:,1), edgee(2,2) == edge(:,2)],2) == 2);
    id3 = find(sum([edgee(3,1) == edge(:,1), edgee(3,2) == edge(:,2)],2) == 2);
    index(e,:) = [t(e,1) midnode(id1) t(e,2) midnode(id2) t(e,3) midnode(id3)];
end
p0 = [p;facecen];
end

