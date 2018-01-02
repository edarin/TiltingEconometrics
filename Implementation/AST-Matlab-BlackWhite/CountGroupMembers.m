function [NG] = CountGroupMembers(G)

G_uniq = unique(G);
N      = length(G_uniq);
NG     = zeros(N,1);

for c = 1:N
    NG(c)   =  length(find(G==G_uniq(c)));
end