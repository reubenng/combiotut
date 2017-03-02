%% Hidden Markov Model for detecting CG rich regions
clear
% phageLambda.fasta genome
pLambda = fileread('phageLambda.fasta');
pLambda = strrep(pLambda,sprintf('\n'),''); % remove \n
symlen = length(pLambda); % observation symbols length

for char = 1:symlen % convert into index for emission [ ATCG ] = [ 1234 ]
    if pLambda(char) == 'A'
        phageLambda(char) = 1;
    elseif pLambda(char) == 'T'
        phageLambda(char) = 2;
    elseif pLambda(char) == 'C'
        phageLambda(char) = 3;
    elseif pLambda(char) == 'G'
        phageLambda(char) = 4;
    end
end

% for testing
% phageLambda = [6 6 6 6 6 6 6 6 6 6 6 6 1 1 1 1 1 1 1 1 2 2 1 2 2 2 2 2 2 2 6 6 6 6 6 6 6 6 6 6 6 6 6]; % for testing
% symlen = length(phageLambda); % observation symbols length

%%
% write all states and probs in the same order
% CG rich region
states = {'L', 'H'}; % hidden states {'AT rich' = (L)ow, 'CG rich' = (H)igh}
numstates = length(states); % number of states
symbols = {'A', 'T', 'C', 'G'}; % number of output symbols for observation
iniprob = [0.5; 0.5]; % initial probabilities of states
iniprob = log(iniprob);
transprob = [0.9998 0.0002; 0.0003 0.9997]; % transition probabilities to another state
transprob = log(transprob);
symprob = [0.2698 0.3237 0.208 0.1985; 0.2459 0.2079 0.2478 0.2984]; % symbol prob given state
symprob = log(symprob);

% % dishonest casino for testing
% states = {'H', 'D'}; % hidden states {'honest', 'dishonest'}
% numstates = length(states); % number of states
% symbols = {'1', '2', '3', '4', '5', '6'}; % number of output symbols for observation
% %         [ H    D ]
% iniprob = [0.5; 0.5]; % initial probabilities of states
% iniprob = log(iniprob);
% %         H [toD stayinH]
% %         D [toH stayinD]
% transprob = [9/10 1/10; 1/10 9/10]; % transition probabilities to another state
% transprob = log(transprob);
% %       H [1    2   3   4   5   6]
% %       D [1    2   3   4   5   6]
% symprob = [1/6 1/6 1/6 1/6 1/6 1/6; 1/10 1/10 1/10 1/10 1/10 1/2]; % symbol prob given state
% symprob = log(symprob);

%%
problist = zeros(numstates, symlen); % state prob for each symbol in sequence
% maxprob is most likely to be in state from max problist
maxprob = zeros(1, symlen); % make a max prob list
likelystate = zeros(numstates, symlen); % most likely state

for i = 1:symlen
    if i == 1 % find initial state
        % ini prob both state * first sym prob both state
        problist(:, i) = iniprob + symprob(:, phageLambda(i)); % prob for both initial state
        % [prob for state, state]
        [~, maxprob(i)] = max(problist(:, i));
        continue;
    end
    
    % starting from i = 2
    for j=1:numstates % for each states
        % next state prob = prob being in each state * trans prob of the current state
        [problist(j, i), likelystate(j, i)] = max(problist(:, i-1) + transprob(j, :)');
    end
    
    % for subsequence symbol after the first
    % state prob = prob of state * prob of sym appearing in state
    problist(:,i) = problist(:, i) + symprob(:, phageLambda(i));
    % predict the state using max prob
    [~, maxprob(i) ]= max(problist(:,i));
end
% maxprob % result without backtracking for testing

%% backward track
[~, maxprob(symlen)] = max(problist(:,symlen)); % finding termination state
statename = cell(1, symlen);
for i=symlen-1:-1:1 % start at second last, back track
    maxprob(i) = likelystate(maxprob(i+1), i+1);
end

for i=1:symlen
    statename{i} = states{maxprob(i)};    
end
statename