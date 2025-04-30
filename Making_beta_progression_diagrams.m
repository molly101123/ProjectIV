% Here we have utilised the following function to create plots of the beta
% model for varying values of beta. Diagrams have been created following
% the method outlined in the examples section of the documentation for the
% function WattsStrogatz, which can be found at:
% https://uk.mathworks.com/help/matlab/math/build-watts-strogatz-small
% -world-graph-model.html#d126e16286.

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


p = linspace(0,1 ,10 ) 
p
%plot(p)
num_p = 10


n = 20
nei = 2
dim = 1
 

g = WattsStrogatz(n,nei,p(10))
plot(g,'NodeColor','k','Layout','circle');
title('Watts-Strogatz Graph with $N = 20$ nodes, $K = 2$, $\beta =1.0000$', ...
    'Interpreter','latex')

p
