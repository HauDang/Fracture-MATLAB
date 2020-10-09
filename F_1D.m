function [ F ] = F_1D( p,bc,bv,Element )
sdof = size(p,1)*2;
F = sparse(sdof,1);
if strcmp(Element,'T6') == 1
    for e = 1:size(bc,1)
        index = bc(e,:);
        dx = max([norm(p(index(3),:) - p(index(1),:)),norm(p(index(3),:) - p(index(2),:)),norm(p(index(2),:) - p(index(1),:))]);
        Fex = [1/6;4/6;1/6]*bv(e,1)*dx;
        Fey = [1/6;4/6;1/6]*bv(e,2)*dx;
        F(index*2-1) = F(index*2-1) + Fex;
        F(index*2) = F(index*2) + Fey;
    end
else
    for e = 1:size(bc,1)
        index = bc(e,:);
        dx = norm(p(index(1),:) - p(index(2),:));
        Fex = [1/2;1/2]*bv(e,1)*dx;
        Fey = [1/2;1/2]*bv(e,2)*dx;
        F(index*2-1) = F(index*2-1) + Fex;
        F(index*2) = F(index*2) + Fey;
    end
end
end
