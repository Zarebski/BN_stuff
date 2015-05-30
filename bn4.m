%-------------------------------------------------------------------------%
% Define the sprinkler BN from bn1.m
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

%-------------------------------------------------------------------------%
% Draw some sample data
%-------------------------------------------------------------------------%

nsamples = 10^3;    % define the number of samples to make
samples = cell(N,nsamples); % preallocate a cell array for the sample data
                            % each row of the array corresponds to an
                            % observation
for i = 1:nsamples
    samples(:,i) = sample_bnet(bnet);
end

data = cell2num(samples);   % since the variables are all scalars we can 
                            % store these samples as a matrix to improve
                            % efficiency

%-------------------------------------------------------------------------%
% Learn the BN parameters
%-------------------------------------------------------------------------%                           
                            
% make a tabula rasa
bnet2 = mk_bnet(dag, node_sizes);   % make a BN with the correct structure
seed = 0;                           % set the seed of the random number                                   
rand('state', seed);                % generator
bnet2.CPD{C} = tabular_CPD(bnet2, C);   % since we don't supply the 
bnet2.CPD{R} = tabular_CPD(bnet2, R);   % parameters tabular_CPD generates
bnet2.CPD{S} = tabular_CPD(bnet2, S);   % random parameters for bnet2
bnet2.CPD{W} = tabular_CPD(bnet2, W);

% estimate the parameters using an increasing number of data points to see
% the learning process
Errors = zeros(1,nsamples);
for j = 1:nsamples
    % find maximum likelihood estimates of the parameters
    bnet3 = learn_params(bnet2, data(:,1:j));

    % veiw the estimated parameters
    CPT3 = cell(1,N);
    for i=1:N
      s3=struct(bnet3.CPD{i});  % violate object privacy
      CPT3{i}=s3.CPT;
      s=struct(bnet.CPD{i});  % violate object privacy
      CPT{i}=s.CPT;
    end

    Errors(j) = norm(CPT{3}-CPT3{3});
end

semilogy(1:nsamples,Errors);
title('Error in parameter MLE');
xlabel('Number of observations used');
ylabel('2-norm of error parameter for "Raining"');
