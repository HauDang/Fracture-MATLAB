function [ id ] = nodesonboundary( p,t,pv,option )
ds = dsegment(p,pv);
% [ds,pv] = point2segment(pv,p);
% error = sum(sum(sqrt((ds-ds0).^2)))
bc = find(min(ds,[],2) <= eps*1e5);
dis = sqrt((p(bc,1) -  pv(1,1)).^2 + (p(bc,2) -  pv(1,2)).^2);
[value,index] = sort(dis);
bcs = bc(index);
if option == 0
    id = bcs;
else
    bar = unique([t(:,[1 2]);t(:,[1 3]);t(:,[2 3])],'rows');
    id1 = bcs(1);id = id1;
    while length(bcs) > 1
        [m,n] = find(bar == id1);
        nodex = bar(m,:);nodex = setdiff(unique(nodex(:)),id1);
        id2 = intersect(nodex,bcs);
        id = [id;id2(1)];
        bcs = setdiff(bcs,id1);
        id1 = id2(1);
    end
end

