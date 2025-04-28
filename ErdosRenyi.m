function A = ErdosRenyi(N,p)
% A = ErdosRenyi(N,p) returns an Erdős–Rényi model random graph
% with N nodes, and p the probability of connection of every possible edge.

% A is adjacency matrix

% To get average degree k, set p = k/(N-1). Exact target not guaranteed due
% to stochasticity.

A = zeros(N,N);
for v = 1:N
    for w = v+1:N
        if rand() < p            
            A(v,w) = 1;
            A(w,v) = 1;
        end
    end
end
end