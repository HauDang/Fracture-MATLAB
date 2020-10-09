function [ p,t ] = RemoveDuplicateNode( p,t )
id = [];
for k = size(p,1):-1:1
    pk = p(k,:);
    local = find(sum([abs(pk(1) - p(1:k-1,1)) <= eps*100 abs(pk(2) - p(1:k-1,2)) <= eps*100],2) == 2);
    if isempty(local) == 0
        id = [id;k local];
    end
end
if isempty(id) == 0
    p(id(:,1),:) = [];
    for n = 1:size(id,1)
        [a,b] = find(t == id(n,1));
        for m = 1:length(a)
            t(a(m),b(m)) = id(n,2);
        end
    end
end
[ t ] = SortingElementIndex( t );
end

