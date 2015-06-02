n = 5; 				% A and B take values 1,...,n
N = 3; 				% the number of nodes
dag = zeros(N, N); 	% preallocat the adj-matrix
A = 1; B = 2; C = 3;
dag([A B], C) = 1; 	% add arcs to out graph

draw_graph(dag, {'A','B','C'}); % produce a figure of the DAG

discrete_nodes = 1:N; 	% list of the discrete nodes
node_sizes = [n n 2]; 	% list of the number of states

% use mk_bnet to make the BN
bnet = mk_bnet(dag, node_sizes, 'discrete', discrete_nodes); 

% put uniform distributions on the variables A and B
bnet.CPD{A} = tabular_CPD(bnet, A, (1/n)*ones(1,n) );
bnet.CPD{B} = tabular_CPD(bnet, B, (1/n)*ones(1,n) );

% C is the identity for {A+B=n} so its CPT is a little more involves
Ccpt = zeros(n,n,2);
for i = 1:n
	for j = 1:n
		if i + j == n
			% constraint is satisfied
			Ccpt(i,j,1) = 0; 	% since 1 is false
			Ccpt(i,j,2) = 1; 	% and 2 is true in bnt
		else
			% constraint is not satisfied 
			Ccpt(i,j,1) = 1; 	
			Ccpt(i,j,2) = 0;
		end
	end
end
bnet.CPD{C} = tabular_CPD(bnet, C, 'CPT', Ccpt);

% choose an inference engine
engine = jtree_inf_engine(bnet); % junction tree inference is exact

% produce a visualisation of the prior joint distribution of A and B 
evidence = cell(1,N);
[engine, ~] = enter_evidence(engine, evidence);
m = marginal_nodes(engine, [A,B]);
fprintf('The prior joint distribution of A and B\n');
disp(m.T);

% instantiate C to reflect constraint being satisfied
evidence{C} = 2; 	
[engine, ~] = enter_evidence(engine, evidence);

% produce a visualisation of the posterior joint distribution of A and B 
m = marginal_nodes(engine, [A B]);
fprintf('The posterior joint distribution of A and B\n');
disp(m.T);
