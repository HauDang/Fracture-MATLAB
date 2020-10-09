function [ IndOut,NodDel,EleDel ] = regionremesh( p1,r,L,p,t,LocalA1,LocalA2,bc )
p00 = p1;
vecx = p00(2:end,1) - p00(1,1); vecy = p00(2:end,2) - p00(1,2);
len = sqrt(vecx.^2 + vecy.^2);
p00(2:end,1) = p00(1,1)+vecx./len*1.5*L;
p00(2:end,2) = p00(1,2)+vecy./len*1.5*L;
L0 = 1.1*min(L,r);
p00(2,:) = p00(1,:) + [vecx(1)/len(1)*L0, vecy(1)/len(1)*L0];
p00(end,:) = p00(1,:) + [vecx(end)/len(end)*L0, vecy(end)/len(end)*L0];
p00(1,:) = p00(end,:);
hold on
plot(p00(:,1), p00(:,2),'r-o')

[in, on] = inpolygon(p(:,1),p(:,2),p00(:,1),p00(:,2)); 
NodClo = union(find(in == 1), find(on == 1));

% LocalO = Oindex
% EleIn = find(any(t == LocalO,2) == 1);
% IndIn = t(EleIn,:);
% IndIn = setdiff(unique(IndIn(:)),LocalO);
tt = t;
% p0 = p(LocalO,:);
% NodClo = find(sqrt((p(:,1) - p0(1)).^2 + (p(:,2) - p0(2)).^2) <= 1.2*L); 
% NodClo = unique(union(NodClo,IndIn));
EleClo = find(any([ismember(t(:,1),NodClo) ismember(t(:,2),NodClo) ismember(t(:,3),NodClo)],2) == 1);
NodAro = t(EleClo,:); NodAro = unique(NodAro(:)); %hold on; plot(p(NodAro,1),p(NodAro,2),'go')
EleDel = find(sum([ismember(t(:,1),NodAro) ismember(t(:,2),NodAro) ismember(t(:,3),NodAro)],2) == 3);
tt(EleDel,:) = []; tt = unique(tt(:)); NodOut = intersect(NodAro,tt); %hold on; plot(p(NodOut,1),p(NodOut,2),'r*')
NodDel = setdiff(NodAro,NodOut); % hold on; plot(p(NodDel,1),p(NodDel,2),'g*')
bar = unique([t(:,[1 2]);t(:,[1 3]);t(:,[2 3])],'rows');NodOut = unique([NodOut;LocalA1;LocalA2]);
NodOut = setdiff(NodOut,[ LocalA1;LocalA2 ]);
IndOut = [LocalA1;LocalA2];
while isempty(NodOut)==0
    [i,j] = find(bar == LocalA1);
    nodex = bar(i,:);nodex = setdiff(nodex(:),LocalA1);
    LocalA1 = intersect(NodOut,nodex);
    if length(LocalA1) >= 2
        id = intersect(LocalA1,bc);
        [i,j] = find(bar == id);
        nod2 = bar(i,:);nod2 = setdiff(unique(nod2(:)),id);
        LocalA1 = intersect(NodOut,nod2);
        NodOut = setdiff(NodOut,[LocalA1;id]);
        IndOut = [IndOut(1:end-1);id;LocalA1;IndOut(end)];
    else    
        IndOut = [IndOut(1:end-1);LocalA1;IndOut(end)];
        NodOut = setdiff(NodOut,LocalA1);
    end
end
end
% hold on
% trisurf(t(eledelete,:),p(:,1),p(:,2),p(:,1)*0,'facecolor','w');