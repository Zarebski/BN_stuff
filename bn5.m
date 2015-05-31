%-------------------------------------------------------------------------%
% As in bn4.m define the sprinkler BN and generate some sample data
%-------------------------------------------------------------------------%

N = 4;                          
dag = zeros(N,N);               
C = 1; S = 2; R = 3; W = 4;     
dag(C,[R S]) = 1;               
dag(R,W) = 1;
dag(S,W)=1;

discrete_nodes = 1:N;           
node_sizes = 2*ones(1,N);       

bnet = mk_bnet(dag, node_sizes);

bnet.CPD{C} = tabular_CPD(bnet, C, [0.5 0.5]);
bnet.CPD{R} = tabular_CPD(bnet, R, [0.8 0.2 0.2 0.8]);
bnet.CPD{S} = tabular_CPD(bnet, S, [0.5 0.9 0.5 0.1]);
bnet.CPD{W} = tabular_CPD(bnet, W, [1 0.1 0.1 0.01 0 0.9 0.9 0.99]);

nsamples = 10^3;    % define the number of samples to make
samples = cell(N,nsamples); % preallocate a cell array for the sample data
                            % each row of the array corresponds to an
                            % observation
for i = 1:nsamples
    samples(:,i) = sample_bnet(bnet);
end

%-------------------------------------------------------------------------%
% Generate some partially observed samples
%-------------------------------------------------------------------------%

samples2 = samples;             % make a copy of the samples
hide = rand(N, nsamples) > 0.5; % generate a matrix of 0's and 1's of the 
                                % same dimensions
[I,J]=find(hide);               % use find to generate coordinate vectors 
                                % for the variables we will ``unobserve"
for k=1:length(I)
  samples2{I(k), J(k)} = [];    % remove the observation of this variable
end
% the braces have been used above since samples2 is a cell array

%-------------------------------------------------------------------------%
% Specify an inference engine and use EM algorithm to learn parameters
%-------------------------------------------------------------------------%

engine2 = jtree_inf_engine(bnet2);
max_iter = 10;
[bnet4, LLtrace] = learn_params_em(engine2, samples2, max_iter);

plot(LLtrace, 'x-');


















