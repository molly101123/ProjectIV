% Here we construct a generalised graph of the clustering and
% characteristic path length of the beta model utilising the function 
% WattsStrogatz as given below. The documentation for this function can be
% found at: 
% https://uk.mathworks.com/help/matlab/math/build-watts-strogatz-small
% -world-graph-model.html#d126e16286

% here we also notably utilise the built in MatLab functions of numnodes,
% neighbours, subgraph, numedges and degree.

% Copyright 2015 The MathWorks, Inc.

function h = WattsStrogatz(N,K,beta)
% H = WattsStrogatz(N,K,beta) returns a Watts-Strogatz model graph with N
% nodes, N*K edges, mean node degree 2*K, and rewiring probability beta.
%
% beta = 0 is a ring lattice, and beta = 1 is a random graph.

% Connect each node to its K next and previous neighbors. This constructs
% indices for a ring lattice.
s = repelem((1:N)',1,K);
t = s + repmat(1:K,N,1);
t = mod(t-1,N)+1;

% Rewire the target node of each edge with probability beta
for source=1:N    
    switchEdge = rand(K, 1) < beta;
    
    newTargets = rand(N, 1); % this strictly differs from the Watts Strogatz model
    newTargets(source) = 0;
    newTargets(s(t==source)) = 0;
    newTargets(t(source, ~switchEdge)) = 0;
    
    [~, ind] = sort(newTargets, 'descend');
    t(source, switchEdge) = ind(1:nnz(switchEdge));
end

h = graph(s,t);
end

function h = clustering(G)
% first I need to identify the total number of nodes in the graph 
local_clustering = repmat([0],1,numnodes(G)) %linspace(0,1,numnodes(G))
for i = 1:numnodes(G)
    neigh = neighbors(G,i)
    sg = subgraph(G,neigh)
    total_edges = numedges(sg)

    deg_vertex = degree(G,i)

    local_clustering(i) = total_edges/nchoosek(deg_vertex,2)
end
total_clustering = mean(local_clustering)
h = total_clustering
end




num_p = 14;
p = linspace(-4,0 ,num_p); 
p = 10.^(p);
iter =20; 
n = 1000;
nei = 5;
dim = 1;
tracking = repmat([0],iter,num_p);
 
for i = 1:num_p 
    for j = 1:iter
        g = WattsStrogatz(n,nei,p(i));
        tracking(j,i) = clustering(g);
    end
end
tracking2 = mean(tracking); % this doesn't work when we only have 1 iteration
max(tracking2)
tracking3 = tracking2/max(tracking2)

y1 = tracking3
x = log10(p)
%plot(log10(p),tracking3,'-o')


num_p = 14;
p = linspace(-4,0 ,num_p); 
p = 10.^(p);
iter =50; %trying to make it more smooth 
n = 1000;
nei = 5;
dim = 1;
char_paths = repmat([0],iter,num_p)
for j = 1:iter
    for i = 1:num_p 
        g = WattsStrogatz(n,nei,p(i))
        char_paths(j,i) = mean(mean(distances(g)))
    end
end
char_paths
char_paths2 = mean(char_paths)
char_paths3 = char_paths2/max(char_paths2)
y2= char_paths3
%plot(log10(p),char_paths3,'-o')


figure
plot(x,y1,'-o')

hold on 

plot(x,y2,'-o')
hold off