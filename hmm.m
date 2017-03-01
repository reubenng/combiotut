%% Hidden Markov Model for detecting CG rich regions
clear
% phageLambda.fasta genome
% phageLambda = fileread('phageLambda.fasta');
% phageLambda = strrep(phageLambda,sprintf('\n'),''); % remove \n

phageLambda = [6 6 6 1 1 1 6 6 6 1]; % for testing

%%
% write all states and probs in the same order
% states = {'X', 'Y'}; % hidden states {'AT rich', 'CG rich'}
% numstates = length(states); % number of states
% symbols = {'A', 'T', 'C', 'G'}; % number of output symbols for observation
% iniprob = [0.5 0.5]'; % initial probabilities of states
% transprob = [0.9998 0.0002; 0.9997 0.0003]; % transition probabilities to another state
% symprob = [0.2698 0.3237 0.208 0.1985; 0.2459 0.2079 0.2478 0.2984]'; % symbol prob given state

% casino for testing
states = {'H', 'D'}; % hidden states {'honest', 'dishonest'}
numstates = length(states); % number of states
symbols = {'1', '2', '3', '4', '5', '6'}; % number of output symbols for observation
% [ H D ]
iniprob = [0.5; 0.5]; % initial probabilities of states
% H [stayinH toD]
% D [stayinD toH]
transprob = [9/10 1/10; 9/10 1/10]; % transition probabilities to another state
% H [1 2 3 4 5 6]
% D [1 2 3 4 5 6]
symprob = [1/6 1/6 1/6 1/6 1/6 1/6; 1/10 1/10 1/10 1/10 1/10 1/2]; % symbol prob given state

symlen = length(phageLambda); % observation symbols length
problist = zeros(numstates, symlen); % all prob for each obs in each state
maxlist = zeros(numstates, symlen); %
maxprob = zeros(1, symlen); % make a max prob list

for i = 1:symlen
    if i == 1 % find initial state
        problist(:, i) = iniprob .* symprob(:, phageLambda(i)); % prob for both initial state
        % [prob for state, state]
        [maxstateprob, maxprob(i)] = max(problist(:, i));
        continue;
    end
    
    for j=1:numstates % for both states
%         [problist(j,i) maxlist(j,i)] = max(problist(:, i-1)' .* transprob(j, :));
        [problist(j, i) maxlist(j, i)] = max(problist(:, i-1) .* transprob(j, :)');
    end
%     problist
%     maxlist
    symprob(:, phageLambda(i))
%     problist(:,i) = problist(:, i) .* symprob(phageLambda(i), :)';
%     [maxstateprob maxprob(i)]=max(problist(:,i));
end

% dec_state=zeros(1,symlen);
% decode_state=cell(1,symlen);
% [pstar dec_state(symlen)]=max(problist(:,symlen));
% 
% for i=symlen-1:-1:1
%     dec_state(i)=maxlist(dec_state(i+1),i+1);
% end
% 
% for i=1:symlen
% decode_state{i}=states{dec_state(i)};    
% end
% decode_state