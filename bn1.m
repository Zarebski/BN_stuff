
%-------------------------------------------------------------------------%
% Create adjacency matrix
%-------------------------------------------------------------------------%

N = 4;                          % number of nodes
dag = zeros(N,N);               % preallocate for the adjacency matrix
C = 1; S = 2; R = 3; W = 4;     % label the nodes
dag(C,[R S]) = 1;               % add arcs into the adjacency matrix
dag(R,W) = 1;
dag(S,W)=1;

%-------------------------------------------------------------------------%
% Visualise the DAG
%-------------------------------------------------------------------------%

% this bnt function allows us to easily visualise this BN
draw_graph(dag, {'Cloudy','Sprinkler','Rain','Wet grass'});
% and put a label on out plot
title('Hello world!');

% the bnt function make_layout returns coordinate vectors for the graph 
% produced by draw_graph  
% [X,Y] = make_layout(dag);

%-------------------------------------------------------------------------%
% Create a BN shell
%-------------------------------------------------------------------------%

discrete_nodes = 1:N;           % decide which nodes will be discrete RVs
node_sizes = 2*ones(1,N);       % decide the number of states they have

% make the BN using the mk_bnet command
% add names to the nodes using the optional 'names' arguement
% use the 'discrete' arguement to specify which variables are discrete
bnet = mk_bnet(dag, node_sizes,'names', {'cloudy','S','R','W'},...
    'discrete', discrete_nodes);

% define the conditional probability tables for each of the variables
bnet.CPD{C} = tabular_CPD(bnet, C, [0.5 0.5]);
bnet.CPD{R} = tabular_CPD(bnet, R, [0.8 0.2 0.2 0.8]);
bnet.CPD{S} = tabular_CPD(bnet, S, [0.5 0.9 0.5 0.1]);
bnet.CPD{W} = tabular_CPD(bnet, W, [1 0.1 0.1 0.01 0 0.9 0.9 0.99]);
% note that the specification of this table requires the user to know how
% the table is constructed STOP for a discussion of the ordering of the
% entries in these vectors see the documentation at
% http://bnt.googlecode.com/svn/trunk/docs/usage.html
% there is also some code there for the assignment of random parameters

%-------------------------------------------------------------------------%
% Inference
%-------------------------------------------------------------------------%

% choose an inference engine
engine = jtree_inf_engine(bnet);    % junction tree inference is exact

% define some evidence
evidence = cell(1,N);       % evidence is that variable 'Wet' is
evidence{W} = 2;            % observed to be in state 2 (true)
% an equivalent potentially more transperant way to say this is 
% evidence{2}=2                           

% feed the evidence into the engine
[engine, ~] = enter_evidence(engine, evidence);

%engine is an updated form of the engine which now incorporates the
%evidence and the ~ is the (unused) log-likelihood of the evidence

% compute the probability that the sprinkler is on given that the grass is
% wet

marg = marginal_nodes(engine, S);   % compute the marginal distribution of
                                    % the variable S 
p_g = marg.T; % define p as the prob the sprinkler is on

fprintf('P(Srinkler=on|Grass=wet) = %f\n',p_g(2));

% compute the probability given the additional evidence that it is raining
evidence{R} = 2;
[engine, loglik] = enter_evidence(engine, evidence);
marg = marginal_nodes(engine, S);
p_gr = marg.T;

fprintf('P(Srinkler=on|Grass=wet,Raining=true) = %f\n',p_gr(2));

% visualise the probability distributions with and evidence on raining
figure;
subplot(1,2,1), bar(p_g);
title('Given grass wet'), axis([0.5,2.5,0,1]);
subplot(1,2,2), bar(p_gr);
title('Given grass wet and raining'), axis([0.5,2.5,0,1]);
close

%-------------------------------------------------------------------------%
% Joint distributions
%-------------------------------------------------------------------------%

evidence = cell(1,N);   % reset evidence
[engine, ~] = enter_evidence(engine, evidence);
m = marginal_nodes(engine, [S R W]);    % marginal distributions of S, R, W

% examples of looking up values in the joint distribution
fprintf('P(Sprinkler=off,Raining=false,Grass=wet) = %f\n',m.T(1,1,2));
fprintf('P(Sprinkler=on,Raining=false,Grass=wet) = %f\n',m.T(2,1,2));

% update the evidence due to observation that Raining=true
evidence{R} = 2;
[engine, ~] = enter_evidence(engine, evidence);
m = marginal_nodes(engine, [S R W], 1); % the 1 needs to be added to get 
                                        % the full table of values

% look at the same examples
fprintf('P(Sprinkler=off,Raining=false,Grass=wet|Raining=true) = %f\n',m.T(1,1,2));
fprintf('P(Sprinkler=on,Raining=false,Grass=wet|Raining=true) = %f\n',m.T(2,1,2));

