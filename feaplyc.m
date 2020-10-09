function[kk,ff] = feaplyc(kk,ff,bcdof,bcvall)

[n,m]=size(bcdof);n = max(n,m);
sdof=size(kk);

for i=1:n
    c=bcdof(i);
    kk(c,:)=0;
    kk(c,c)=1;
    ff(c,1)=bcvall(i);
end
end