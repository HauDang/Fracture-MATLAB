function [ K ] = K_FEM( D,p,t,CQPE,Element )
sdof = size(p,1)*2;
K = sparse(sdof,sdof);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (strcmp(Element,'T3') == 1 )
    for e = 1:size(t,1)
        X = p(t(e,:),1);
        Y = p(t(e,:),2);
        index = [t(e,1)*2-1, t(e,1)*2, t(e,2)*2-1, t(e,2)*2, t(e,3)*2-1, t(e,3)*2];
        [ dNdx, dNdy, Ae, detJac ] = T3_Element ( X,Y );
        B = [dNdx(1) 0,dNdx(2) 0,dNdx(3) 0;...
             0 dNdy(1),0 dNdy(2),0 dNdy(3);...
             dNdy(1) dNdx(1),dNdy(2) dNdx(2),dNdy(3) dNdx(3)];
        Ke = B'*D*B*Ae;
        K(index,index) = K(index,index) + Ke;
    end
elseif (strcmp(Element,'T6') == 1 )
    for e = 1:size(t,1)
        X = p(t(e,:),1);
        Y = p(t(e,:),2);
        if (X(2) - X(1))*(Y(end) - Y(1)) - (Y(2) - Y(1))*(X(end) - X(1)) < 0
            t(e,:) = t(e,[1 end:-1:2]);
            X = p(t(e,:),1);
            Y = p(t(e,:),2);
        end
        index = [t(e,1)*2-1, t(e,1)*2, t(e,2)*2-1, t(e,2)*2, t(e,3)*2-1, t(e,3)*2, t(e,4)*2-1, t(e,4)*2, t(e,5)*2-1, t(e,5)*2, t(e,6)*2-1, t(e,6)*2];
        Ke = 0;
        if ismember(e,CQPE) == 0
            [poi,wei] = GaussPointRule(6,'T3');
            [ N, dNdx, dNdy, detJac ] = T6_Element( X,Y,poi(:,1),poi(:,2) );
            for i = 1:size(poi,1)
                B = [dNdx(i,1) 0,dNdx(i,2) 0,dNdx(i,3) 0,dNdx(i,4) 0,dNdx(i,5) 0,dNdx(i,6) 0;...
                    0 dNdy(i,1),0 dNdy(i,2),0 dNdy(i,3),0 dNdy(i,4),0 dNdy(i,5),0 dNdy(i,6);...
                    dNdy(i,1) dNdx(i,1),dNdy(i,2) dNdx(i,2),dNdy(i,3) dNdx(i,3),dNdy(i,4) dNdx(i,4),dNdy(i,5) dNdx(i,5),dNdy(i,6) dNdx(i,6)];
                Ke = Ke + B'*D*B*detJac(i)*wei(i);
            end
        else
            [poi,wei] = GaussPointRule(6,'Q4');
            for i = 1:size(poi,1)
                xis = poi(i,1);
                eta = poi(i,2);
                [ dNdx, dNdy, detJac ] = Q8_Element( X,Y,xis,eta );
                B = [dNdx(1)+dNdx(7)+dNdx(8) 0,dNdx(2) 0,dNdx(3) 0,dNdx(4) 0,dNdx(5) 0,dNdx(6) 0;...
                    0 dNdy(1)+dNdy(7)+dNdy(8),0 dNdy(2),0 dNdy(3),0 dNdy(4),0 dNdy(5),0 dNdy(6);...
                    dNdy(1)+dNdy(7)+dNdy(8) dNdx(1)+dNdx(7)+dNdx(8),dNdy(2) dNdx(2),dNdy(3) dNdx(3),dNdy(4) dNdx(4),dNdy(5) dNdx(5),dNdy(6) dNdx(6)];
                Ke = Ke + B'*D*B*detJac*wei(i);
            end
        end
        K(index,index) = K(index,index) + Ke;
    end
end
end


function [ N, dNdx, dNdy, detJac ] = T6_Element( X,Y,xis,eta )
N = [(1 - xis - eta).*(1 - 2*xis - 2*eta),...
    4*xis.*(1 - xis - eta),...
    xis.*(2*xis - 1),...
    4*xis.*eta,...
    eta.*(2*eta - 1),...
    4*eta.*(1 - xis - eta)];
% x = N(:,1)*X(1) + N(:,2)*X(2) + N(:,3)*X(3) + N(:,4)*X(4) + N(:,5)*X(5) + N(:,6)*X(6);
% y = N(:,1)*Y(1) + N(:,2)*Y(2) + N(:,3)*Y(3) + N(:,4)*Y(4) + N(:,5)*Y(5) + N(:,6)*Y(6);
dNdxis = [4*eta + 4*xis - 3,...
    4 - 8*xis - 4*eta,...
    4*xis - 1,...
    4*eta,...
    0*xis,...
    -4*eta];
dNdeta = [4*eta + 4*xis - 3,...
    -4*xis,...
    0*xis,...
    4*xis,...
    4*eta - 1,...
    4 - 4*xis - 8*eta];

dxdxis = dNdxis*X;
dxdeta = dNdeta*X;
dydxis = dNdxis*Y;
dydeta = dNdeta*Y;

detJac = dxdxis.*dydeta-dxdeta.*dydxis;
dxisdx =  dydeta./detJac;
dxisdy = -dxdeta./detJac;
detadx = -dydxis./detJac;
detady =  dxdxis./detJac;
for i = 1:length(dxisdx)
    dNdx(i,:) = dNdxis(i,:)*dxisdx(i) + dNdeta(i,:)*detadx(i);
    dNdy(i,:) = dNdxis(i,:)*dxisdy(i) + dNdeta(i,:)*detady(i);
end
end

function [ dNdx, dNdy, Ae, detJac ] = T3_Element ( X,Y )
x1=X(1);x2=X(2);x3=X(3);
y1=Y(1);y2=Y(2);y3=Y(3);
dN1dx= (y2 - y3)/(x3*y1 - x2*y1 + x2*y3 - x3*y2 + x1*(y2 - y3));
dN1dy=-(x2 - x3)/(x3*y1 - x2*y1 + x2*y3 - x3*y2 + x1*(y2 - y3));
dN2dx= (y1 - y3)/(x1*y3 - x1*y2 - x3*y1 + x3*y2 + x2*(y1 - y3));
dN2dy= (x1 - x3)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2);
dN3dx= (y1 - y2)/(x3*y1 - x2*y1 + x2*y3 - x3*y2 + x1*(y2 - y3));
dN3dy=-(x1 - x2)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2);
l1=norm([x2-x1;y2-y1]);
l2=norm([x3-x1;y3-y1]);
l3=norm([x2-x3;y2-y3]);
s=(l1+l2+l3)/2;
Ae=sqrt(s*(s-l1)*(s-l2)*(s-l3));
detJac = Ae*2;
dNdx=[dN1dx,dN2dx,dN3dx];
dNdy=[dN1dy,dN2dy,dN3dy];
end

function [ dNdx, dNdy, detJac ] = Q8_Element( X,Y,xis,eta )
poi = [-1 -1;...
    0 -1;...
    1 -1;...
    1 0;...
    1 1;...
    0 1;...
    -1 1;...
    -1 0];
for i = 1:8
    xisi = poi(i,1); etai = poi(i,2);
    N(i) = ( (1 + xis*xisi)*(1 + eta*etai) - (1 - xis^2)*(1 + eta*etai) - (1 - eta^2)*(1 + xis*xisi) )*xisi^2*etai^2/4 ...
        + (1 - xis^2)*(1 + eta*etai)*(1 - xisi^2)*etai^2/2 + (1 - eta^2)*(1 + xis*xisi)*(1 - etai^2)*xisi^2/2;
    dNdxis(i) = ( xisi*(1 + eta*etai) + 2*xis*(1 + eta*etai) - xisi*(1 - eta^2) )*xisi^2*etai^2/4 ...
        - 2*xis*(1 + eta*etai)*(1 - xisi^2)*etai^2/2 + xisi*(1 - eta^2)*(1 - etai^2)*xisi^2/2;
    dNdeta(i) = ( etai*(1 + xis*xisi) - etai*(1 - xis^2) + 2*eta*(1 + xis*xisi) )*xisi^2*etai^2/4 ...
        + etai*(1 - xis^2)*(1 - xisi^2)*etai^2/2 - 2*eta*(1 + xis*xisi)*(1 - etai^2)*xisi^2/2;
end
% x = (N(1) + N(7) + N(8))*X(1) + N(2)*X(2) + N(3)*X(3) + N(4)*X(4) + N(5)*X(5) + N(6)*X(6);
% y = (N(1) + N(7) + N(8))*Y(1) + N(2)*Y(2) + N(3)*Y(3) + N(4)*Y(4) + N(5)*Y(5) + N(6)*Y(6);
dxdxis = (dNdxis(1) + dNdxis(7) + dNdxis(8))*X(1) + dNdxis(2)*X(2) + dNdxis(3)*X(3) + dNdxis(4)*X(4) + dNdxis(5)*X(5) + dNdxis(6)*X(6);
dxdeta = (dNdeta(1) + dNdeta(7) + dNdeta(8))*X(1) + dNdeta(2)*X(2) + dNdeta(3)*X(3) + dNdeta(4)*X(4) + dNdeta(5)*X(5) + dNdeta(6)*X(6);
dydxis = (dNdxis(1) + dNdxis(7) + dNdxis(8))*Y(1) + dNdxis(2)*Y(2) + dNdxis(3)*Y(3) + dNdxis(4)*Y(4) + dNdxis(5)*Y(5) + dNdxis(6)*Y(6);
dydeta = (dNdeta(1) + dNdeta(7) + dNdeta(8))*Y(1) + dNdeta(2)*Y(2) + dNdeta(3)*Y(3) + dNdeta(4)*Y(4) + dNdeta(5)*Y(5) + dNdeta(6)*Y(6);
detJac = dxdxis*dydeta-dxdeta*dydxis;
dxisdx =  dydeta/detJac;
dxisdy = -dxdeta/detJac;
detadx = -dydxis/detJac;
detady =  dxdxis/detJac;
dNdx = dNdxis*dxisdx + dNdeta*detadx;
dNdy = dNdxis*dxisdy + dNdeta*detady;
end

