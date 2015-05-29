
%-------------------------------------------------------------------------%
% First create an adjacency matrix
%-------------------------------------------------------------------------%

N = 4;                          % number of nodes
dag = zeros(N,N);               % preallocate for the adjacency matrix
C = 1; S = 2; R = 3; W = 4;     % label the nodes
dag(C,[R S]) = 1;               % add arcs into the adjacency matrix
dag(R,W) = 1;
dag(S,W)=1;

%-------------------------------------------------------------------------%
% Now we can visualise this graph
%-------------------------------------------------------------------------%

% this bnt function allows us to easily visualise this BN
draw_graph(dag, {'Cloudy','Sprinkler','Rain','Wet grass'});
% and put a label on out plot
title('Hello world!');

% the bnt function make_layout returns coordinate vectors for the graph 
% produced by draw_graph  
% [X,Y] = make_layout(dag);


