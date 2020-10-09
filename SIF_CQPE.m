function [ Keq,thetha,L ] = SIF_CQPE( E,nuy,k,U,p,t,TipEle,Element )
% TipEle: Elements around Tip
G = E/(2*(1 + nuy));
if size(t,2) == 6
    e = t(TipEle(1),3);
    d = t(TipEle(1),2);
    a = t(TipEle(1),1);
    b = t(TipEle(end),6);
    c = t(TipEle(end),5);
else
    e = t(TipEle(1),2);
    a = t(TipEle(1),1);
    c = t(TipEle(end),3);
end
L1 = norm([p(c,1) - p(a,1),p(c,2) - p(a,2)]);
L2 = norm([p(e,1) - p(a,1),p(c,2) - p(a,2)]);
L = 1/2*(L1 + L2);

B1 = p(e,:);
B2 = p(c,:);
O = p(a,:);
B = 1/2*(B1+B2);
[ angle ] = anglebetween3points( [1 0],[0 0],O-B );
x1 = 1; y1 = 0;
x2 = 0; y2 = 0;
x3 = O(1) - B(1); y3 = O(2) - B(2);
d0 = (x2 - x1)*(y3 - y1) - (y2 - y1)*(x3 - x1);
if abs(angle) > eps && abs(abs(angle) - pi) > eps
    angle0 = 2*pi - sign(d0)*angle;
else
    angle0 = angle;
end

Trans = [cos(angle0) sin(angle0);-sin(angle0) cos(angle0)];
Ux = U(1:2:end);
Uy = U(2:2:end);

if strcmp(Element,'T6')
    pc = Trans*[Ux(c);Uy(c)];
    pb = Trans*[Ux(b);Uy(b)];
    pe = Trans*[Ux(e);Uy(e)];
    pd = Trans*[Ux(d);Uy(d)];
    K1 = E/3/(1+k)/(1+nuy)*sqrt(2*pi/L)*(4*(pb(2) - pd(2)) - 0.5*(pc(2) - pe(2)));
    K2 = E/3/(1+k)/(1+nuy)*sqrt(2*pi/L)*(4*(pb(1) - pd(1)) - 0.5*(pc(1) - pe(1)));
elseif strcmp(Element,'T3')
    pc = Trans*[Ux(c);Uy(c)];
    pe = Trans*[Ux(e);Uy(e)];
    K1 = E/2/(1+k)/(1+nuy)*sqrt(2*pi/L)*(pc(2) - pe(2));
    K2 = E/2/(1+k)/(1+nuy)*sqrt(2*pi/L)*(pc(1) - pe(1));
end

alpha = K1/K2;
ang = [2*atand(1/4*(alpha + sqrt(alpha^2 + 8)));2*atand(1/4*(alpha - sqrt(alpha^2 + 8)))];
d2S = -K1*(cosd(ang/2) + 3*cosd(3*ang/2)) + K2*(sind(ang/2) + 9*sind(3*ang/2));
id = find(d2S<0);
if isempty(id)
    thetha = 0;
else
%     thetha = ang(id);
    thetha = sign(ang(id))*min(80,abs(ang(id)));
end
Keq = [K1, K2];


