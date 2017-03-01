%% Hidden Markov Model for detecting CG rich regions
clear
% phageLambda.fasta genome
% phageLambda = fileread('phageLambda.fasta');
% phageLambda = strrep(phageLambda,sprintf('\n'),''); % remove \n

phageLambda = [6 6 6 1 1 1 6 6 6 1]; % for testing

%%
% write all states and probs in the same order
% CG rich region
% states = {'X', 'Y'}; % hidden states {'AT rich', 'CG rich'}
% numstates = length(states); % number of states
% symbols = {'A', 'T', 'C', 'G'}; % number of output symbols for observation
% iniprob = [0.5 0.5]'; % initial probabilities of states
% transprob = [0.9998 0.0002; 0.9997 0.0003]; % transition probabilities to another state
% symprob = [0.2698 0.3237 0.208 0.1985; 0.2459 0.2079 0.2478 0.2984]'; % symbol prob given state

% dishonest casino for testing
states = {'H', 'D'}; % hidden states {'honest', 'dishonest'}
numstates = length(states); % number of states
symbols = {'1', '2', '3', '4', '5', '6'}; % number of output symbols for observation
%         [ H    D ]
iniprob = [0.5; 0.5]; % initial probabilities of states
%         H [toD stayinH]
%         D [toH stayinD]
transprob = [1/10 9/10; 1/10 9/10]; % transition probabilities to another state
%       H [1    2   3   4   5   6]
%       D [1    2   3   4   5   6]
symprob = [1/6 1/6 1/6 1/6 1/6 1/6; 1/10 1/10 1/10 1/10 1/10 1/2]; % symbol prob given state

%%
symlen = length(phageLambda); % observation symbols length
problist = zeros(numstates, symlen); % state prob for each symbol in sequence
% maxprob is most likely to be in state from max problist
maxprob = zeros(1, symlen); % make a max prob list
likelystate = zeros(numstates, symlen); % most likely state

for i = 1:symlen
    if i == 1 % find initial state
        % ini prob both state * first sym prob both state
        problist(:, i) = iniprob .* symprob(:, phageLambda(i)); % prob for both initial state
        % [prob for state, state]
        [~, maxprob(i)] = max(problist(:, i));
        continue;
    end
%     problist
    % starting from i = 2
    for j=1:numstates % for each states
        % next state prob = prob being in each state * trans prob of the current state
%         problist(:, i-1)
        [problist(j, i), likelystate(j, i)] = max(problist(:, i-1) .* transprob(j, :)');
    end
    % likelystate will all be state 2 since it is most probable??
    
    % for subsequence symbol after the first
    % state prob = prob of state * prob of sym appearing in state
    problist(:,i) = problist(:, i) .* symprob(:, phageLambda(i));
    % predict the state using max prob
    [~, maxprob(i)]=max(problist(:,i));
end
% maxprob % result without backtracking for testing

%% backward track
statename=cell(1,symlen); % for storing actual state name
% finding termination state
[~, maxprob(symlen)] = max(problist(:,symlen));

for i=symlen-1:-1:1 % start at second last, back track
    maxprob(i) = likelystate(maxprob(i+1), i+1);
end
% maxprob is base on likelystate and since likelystate is all 2 this will be all state 2 as well

for i=1:symlen
    statename{i} = states{maxprob(i)};    
end
% statename
