%% Smith-Waterman Algorithm
clear
% blosum50 scoring matrix
blosum50mat = dlmread('blosum50.txt');
% blosum50mat
blosum50 = 'ARNDCQEGHILKMFPSTWYV'; % blosum50 matching sequence
cost = -8; % insertion or deletion cost

% Run on HEAGAWGHEE versus PAWHEAE
% sequence1 = 'HEAGAWGHEE'; 
% sequence2 = 'PAWHEAE';

% Compare to lecture
% Match the protein sequence PQPTTPVSSFTSGSMLGRTDTALTNTYSAL with PSPTMEAVEASTASHPHSTSSYFATTYYHL
sequence1 = 'MQNSHSGVNQLGGVFVNGRPLPDSTRQKIVELAHSGARPCDISRILQVSNGCVSKILGRY';
sequence2 = 'TDDECHSGVNQLGGVFVGGRPLPDSTRQKIVELAHSGARPCDISRI';

%% Initialise matrix
s1 = length(sequence1); % find sequence length for making matrix
s2 = length(sequence2);

mat = zeros(s2 + 1, s1 + 1); % fill in all zero
back = cell(s2 + 1, s1 + 1); % for back tracing directions

sx = 0;
sy = 0;
% initialise multiples of -8 score along both axis
for x = 2:s1 + 1 
    back{1, x} = '0'; % value from left
end
for y = 2:s2 + 1
    back{y, 1} = '0'; % value from top
end

%% Fill in matrix
for x = 2:s1 + 1
    for y = 2:s2 + 1
        % find substitution score
        if x <= s1+1 % if index not exceed length of sequence
            xl = sequence1(x-1); % find letters in the position
        end
        xpos = strfind(blosum50, xl); % use letter find pos in blosum50 seq
        % repeat along y-axis
        if y <= s2+1
            yl = sequence2(y-1);
        end
        ypos = strfind(blosum50, yl);
        subscore = blosum50mat(xpos, ypos); % use pos find subscore in blosum mat
        
        score = [mat(y, x-1) + cost, mat(y-1, x) + cost, mat(y-1, x-1) + subscore, 0];
        [mat(y, x), maxcase] = max(score);
        % if either max score same as score from diaganal
        if (mat(y,x-1)+cost)==(mat(y-1,x-1)+subscore) || (mat(y-1,x)+cost)==(mat(y-1,x-1)+subscore)
            back{y, x} = 'd'; % value from diagonal
        else
            switch maxcase
                case 1
                    back{y, x} = 'l'; % value from left
                case 2
                    back{y, x} = 't'; % value from top
                case 3
                    back{y, x} = 'd'; % value from diagonal
                case 4
                    back{y, x} = '0';
            end
        end
    end
end
mat

%% Traceback route
[maxscore, position] = max(mat(:));
[i, j] = ind2sub(size(mat),position);

loop = 1;
l = 1;
% i = s2 + 1;
% j = s1 + 1;
while loop == 1
    if j == 1 && i == 1 % end while loop when reaches 0
        loop = 0;
    end
    if back{i, j} == 'l'
        c = sequence1(j-1);
        match1(l) = c;
        match2(l) = '_';
        j = j - 1;
    elseif back{i, j} == 't'
        match1(l) = '_';
        c = sequence2(i-1);
        match2(l) = c;
        i = i - 1;
    elseif back{i, j} == 'd'
        c = sequence2(i-1); % find letter at that pos
        match2(l) = c; % assign to matching seq
        c = sequence1(j-1);
        match1(l) = c;
        i = i - 1;
        j = j - 1;
    elseif back{i, j} == '0'
        loop = 0;
    end
    l = l + 1; % increment match seq index
end
match1 = fliplr(match1);
match2 = fliplr(match2);
match1
match2

[~, alignment] = swalign(sequence1, sequence2, 'scoringmatrix', 'BLOSUM50')