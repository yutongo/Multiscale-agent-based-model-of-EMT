function PosNormalRandom(m,n,a,b)
# b is standard deviation, a is mean
rv = b.*rand(m,n) + a
ind = find(rv->rv<=0,rv)
while ~isempty(ind)
    for i = 1:length(ind)
        rv[ind[i]] = b.*rand(1) + a;
    end
    ind = find(rv->rv<=0,rv);
end
return rv
end
