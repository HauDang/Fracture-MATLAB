function [ ] = drawmodel( p,t,node,option )
% patch('vertices',p,'faces',t,'edgecol','b','facecol','n','MarkerSize',28,'Marker','.');
patch('vertices',p,'faces',t,'edgecol','b','facecol','n','MarkerSize',6,'Marker','.');
% drawnow;
axis image;axis off
% axis([-0.1 1.1 -0.1 2.2])
%  hold on; plot(p([143 144 145 147 155 146 71 70 69],1),p([143 144 145 147 155 146 71 70 69],2))
 
%  hold on; plot(p([141 142 143 145 153 144 70 69 68],1),p([141 142 143 145 153 144 70 69 68],2))
switch option
    case 1
        hold on
        for i=1:size(p,1)
            h=text(p(i,1),p(i,2),0,int2str(i));
            set(h,'fontsize',16,'color','K');
        end
        axis off
        hold on
        for i=1:size(t,1)
            xx=p(t(i,:),1);
            yy=p(t(i,:),2);
            cenx=sum(xx)/node; ceny=sum(yy)/node;
            h=text(cenx,ceny,0,int2str(i));
            set(h,'fontsize',10,'color','b');
        end
end
end

