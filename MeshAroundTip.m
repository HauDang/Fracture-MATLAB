function [ p,t,B1,B2,Paround,indexCQPE ] = MeshAroundTip( p,t,tip,pv0,L,d0 )
Oindex = find(sum([tip(1) == p(:,1) tip(2) == p(:,2)],2) == 2);
[ bc ] = nodesonboundary( p,t,pv0,0 );
EleIn = find(any(t == Oindex,2) == 1);
IndIn = t(EleIn,:);
IndIn = setdiff(unique(IndIn(:)),Oindex);
[ B1index,B2index ] = findclosetpoints( p,IndIn,Oindex,bc );
%figure; plot(coord(IndOut,1),coord(IndOut,2),'-')
B1 = p(B1index,:); B2 = p(B2index,:); O = p(Oindex,:);
l1 = norm(O-B1);
l2 = norm(O-B2);
if l1 > 1.7*l2
    pnew = 1/2*(O + B1);
    idnew = size(p,1)+1;
    id = [Oindex B1index];
    EleIn = find(sum([any(t == Oindex,2),any(t == B1index,2)],2) == 2);
    id0 = setdiff(t(EleIn,:),id);
    indexnew = [id0 Oindex idnew;id0 B1index idnew];
    p = [p;pnew];
    t(EleIn,:) = [];
    t = [t;indexnew];
    B1index = idnew;
    B1 = p(B1index,:);
elseif l2 > 1.7*l1
    pnew = 1/2*(O + B2);
    idnew = size(p,1)+1;
    id = [Oindex B2index];
    EleIn = find(sum([any(t == Oindex,2),any(t == B2index,2)],2) == 2);
    id0 = setdiff(t(EleIn,:),id);
    indexnew = [id0 Oindex idnew;id0 B2index idnew];
    p = [p;pnew];
    t(EleIn,:) = [];
    t = [t;indexnew];
    B2index = idnew;
    B2 = p(B2index,:);
end

%% quarter point elements: region 1
[p1,t1,inOut,InCen, r] = mesharoundcracktip(B1,B2,O,L);
B1 = p1(inOut(end),:); B2 = p1(inOut(1),:);%hold on;plot(B1(1,1),B1(1,2),'bo')
%% around quater point elemnts: region 2
nodesaroundB1 = t(find(any(t == B1index,2) == 1),:); nodesaroundB1 = unique(nodesaroundB1); A1index = setdiff(intersect(nodesaroundB1,bc),[B1index,Oindex]);
nodesaroundB2 = t(find(any(t == B2index,2) == 1),:); nodesaroundB2 = unique(nodesaroundB2); A2index = setdiff(intersect(nodesaroundB2,bc),[B2index,Oindex]);
[ IndOut,NodDel,EleDel ] = regionremesh( p1,r,L,p,t,A1index,A2index,bc );

pv = [p1(inOut(end),1) p1(inOut(end),2);p(IndOut,1) p(IndOut,2);p1(inOut,1) p1(inOut,2)]; % figure; hold on;plot(pv(:,1),pv(:,2),'-*r')
%% total mesh: region 2
[ p2,t2 ] = GenerationMeshNoExtraNode( pv,d0 );
%% connect region 1 and region 2
[ p12,t12 ] = FixMesh( p2,t2,p1,t1 );

p(NodDel,:) = []; p0 = p(:,[1 2]);
t(EleDel,:) = []; t0 = t;
for k = length(NodDel):-1:1;
    t0(t0>NodDel(k)) = t0(t0>NodDel(k)) - 1;
end
%% connect region 1 2 and region 3
[ p,t ] = FixMesh( p0,t0,p12,t12 );
Paround = p1(2:end,:); 
id0 = find(sum([tip(1) == p(:,1) tip(2) == p(:,2)],2) == 2);
idb = find(sum([B1(1) == p(:,1) B1(2) == p(:,2)],2) == 2);
idt = find(sum([B2(1) == p(:,1) B2(2) == p(:,2)],2) == 2);
indexCQPE = [idb id0 idt];
end
% hold on
% trisurf(t1,p1(:,1),p1(:,2),p1(:,1)*0,'facecolor','g');figure; hold on; drawmodel( p12,t12 )