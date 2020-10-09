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