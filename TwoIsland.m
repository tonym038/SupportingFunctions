function A = TwoIsland(n)
% A = TwoIsland(n) returns a two WS-island graph with n nodes

% A is adjacency matrix

A11 = WS(n/2,10,0.1);
A22 = WS(n/2,10,0.1);
A = blkdiag([A11, zeros(n/2);zeros(n/2), A22]);
A = eye(n)+A;

np=10;
C1 = 1:n/2;
C2 = n/2+1:n;
for t = 1:np
    i = randsample(C1,1); 
    j = randsample(C2,1);
    A(i,j) = 1; A(j,i) = 1;
end   
end