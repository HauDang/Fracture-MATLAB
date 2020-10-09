function [ bc,bv ] = NeuExtension( p,t,parou,tor,EleTyp )
[ nodetop ] = nodesonboundary( p,t,[parou(4,:);parou(3,:)],0 );
[ nodebot ] = nodesonboundary( p,t,[parou(1,:);parou(2,:)],0 );
if strcmp(EleTyp,'T6') == 1
    bct = zeros((length(nodetop)-1)/2,3);
    bvt = zeros((length(nodetop)-1)/2,2);
    
    for i = 1:(length(nodetop)-1)/2
        index = [2*i-1 2*i 2*i+1];
        bct(i,:) = nodetop(index);
        bvt(i,:) = [0, tor];
    end
    
    bcb = zeros((length(nodebot)-1)/2,3);
    bvb = zeros((length(nodebot)-1)/2,2);
    for i = 1:(length(nodebot)-1)/2
        index = [2*i-1 2*i 2*i+1];
        bcb(i,:) = nodebot(index);
        bvb(i,:) = [0, -tor];
    end
else
    bct = zeros(length(nodetop)-1,2);
    bvt = zeros(length(nodetop)-1,2);
    for i = 1:length(nodetop)-1
        index = [i i+1];
        bct(i,:) = nodetop(index);
        bvt(i,:) = [0, tor];
    end
    
    bcb = zeros(length(nodebot)-1,2);
    bvb = zeros(length(nodebot)-1,2);
    for i = 1:length(nodebot)-1
        index = [i i+1];
        bcb(i,:) = nodebot(index);
        bvb(i,:) = [0, -tor];
    end
end
bc = [bct;bcb];
bv = [bvt;bvb];
end

