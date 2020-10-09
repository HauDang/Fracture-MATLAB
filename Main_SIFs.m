clc
% close all
clear
format long
CQPEs = 1;
EleTyp = 'T6';
E = 1; nuy = 0.3;
% plane stress / very thin plate
thick = 0.2;
D = thick*E/(1-nuy^2)*[1 nuy 0;
    nuy 1 0;...
    0 0 (1-nuy)/2];
k = (3 - nuy)/(1 + nuy);
% plane strain / very long beam
% thick = 20;
% D = thick*E*(1 - nuy)/((1+nuy)*(1 - 2*nuy))*[1 nuy/(1 - nuy) 0;
%     nuy/(1 - nuy) 1 0;...
%     0 0 (1-2*nuy)/(2*(1-nuy))];
% k = 3 - 4*nuy;
p3 = load('coordinate.txt'); % node coordinates
t3 = load('element.txt'); % element indices 
fra1 = load('fra1.txt'); % coordinate of fracture faces
fra2 = load('fra2.txt'); % coordinate of fracture faces
tip = {[fra1(1,1), fra1(1,2)],...
       [fra1(end,1), fra1(end,2)]}; % fracture tip
drawmodel( p3,t3,3,0 ); axis image

cbl = intersect(find(p3(:,1) == min(p3(:,1))), find(p3(:,2) == min(p3(:,2))));
cbr = intersect(find(p3(:,1) == max(p3(:,1))), find(p3(:,2) == min(p3(:,2))));
ctl = intersect(find(p3(:,1) == min(p3(:,1))), find(p3(:,2) == max(p3(:,2))));
ctr = intersect(find(p3(:,1) == max(p3(:,1))), find(p3(:,2) == max(p3(:,2))));
parou = p3([cbl cbr ctr ctl cbl],:);
if strcmp(EleTyp,'T6')
    [ p6,t6 ] = T6_pretreatment( p3,t3 );
    CQPE = [];
    for i = 1:length(tip)
        tipi = tip{i};
        Tipindex = find(sum([abs(tipi(1,1) - p3(:,1))<1e5*eps abs(tipi(1,2) - p3(:,2))<1e5*eps],2) == 2);
        qpei = find(any(t3 == Tipindex,2) == 1);
        CQPEi{i} = qpei;
        p6(t6(qpei,2),1) = 3/4*p6(t6(qpei,1),1) + 1/4*p6(t6(qpei,3),1);
        p6(t6(qpei,2),2) = 3/4*p6(t6(qpei,1),2) + 1/4*p6(t6(qpei,3),2);
        p6(t6(qpei,6),1) = 3/4*p6(t6(qpei,1),1) + 1/4*p6(t6(qpei,5),1);
        p6(t6(qpei,6),2) = 3/4*p6(t6(qpei,1),2) + 1/4*p6(t6(qpei,5),2);
        CQPE = [CQPE;qpei];
    end
else
    CQPE = [];
    for i = 1:2
        tipi = tip{i};
        Tipindex = find(sum([abs(tipi(1,1) - p3(:,1))<1e5*eps abs(tipi(1,2) - p3(:,2))<1e5*eps],2) == 2);
        qpei = find(any(t3 == Tipindex,2) == 1);
        CQPEi{i} = qpei;
    end
    p6 = p3; t6 = t3;
end

[ bc,bv ] = NeuExtension( p6,t6,parou,1e-4,EleTyp );
K = K_FEM( D,p6,t6,CQPE,EleTyp );
F = thick*F_1D( p6,bc,bv,EleTyp );
[ nodebot ] = nodesonboundary( p6,t6,[parou(1,:);parou(2,:)],0 );
[ nodetop ] = nodesonboundary( p6,t6,[parou(4,:);parou(3,:)],0 );
udof = [nodebot*2-1; nodebot*2];
uval = udof*0;

[K,F] = feaplyc2(K,F,udof,uval);
U = K\F;
KeqFem = [];
for i = 1:2
    qpei = CQPEi{1,i};
    [ Keq,thetha,L ] = SIF_CQPE( E,nuy,k,U,p6,t6,qpei,EleTyp );
    CraAng(i) = thetha;
    KeqFem = [KeqFem abs(Keq)];
end
KeqFem

bc1 = nodesonboundary( p6,t6,fra1,0 );
bc2 = nodesonboundary( p6,t6,fra2,0 );

ux = U(1:2:size(p6,1)*2);
uy = U(2:2:size(p6,1)*2);

xx = p6(bc2,1);
yy = p6(bc2,1);
yy1 = full(uy(bc1));
yy2 = full(uy(bc2));
plot(xx,yy1,'ro')
hold on
plot(xx,yy2,'bo')
figure
trisurf(t6,p6(:,1) + 1e3*ux,p6(:,2) + 1e3*uy,ux*0);view(0,90); axis image
