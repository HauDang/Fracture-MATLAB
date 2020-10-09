function [ Ae,areamesh ] = computearea( X,Y )
x1=X(:,1);x2=X(:,2);x3=X(:,3);
y1=Y(:,1);y2=Y(:,2);y3=Y(:,3);
l1=sqrt((x2-x1).^2 + (y2-y1).^2);
l2=sqrt((x3-x1).^2 + (y3-y1).^2);
l3=sqrt((x2-x3).^2 + (y2-y3).^2);
s=(l1+l2+l3)./2;
Ae=sqrt(s.*(s-l1).*(s-l2).*(s-l3));
areamesh = sum(Ae);
end

