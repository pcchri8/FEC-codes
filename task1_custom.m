%Team #8
%Christof Petros 9928
%Grigoriou Stergios 9564
%Zacharioudaki Danai 9418

%% task_1 for the requested inputs

% Returns bit error rate, block error rate and number of bits used for each
% combination of M,d and modul_name parameter in a 4-D matrix (M*d*m_n*err) 
% It also returns the data rate given a symbol time Ts (Default is 1) for
% each combination of M and d in a length(M)*length(d) array.

%-IN 
% SNR  scalar
% M,d  vectors
% modul_name cell_array (eg.{'PAM','QAM','PSK'})
%Ts scalar
%-OUT

function [err,data_rates] = task1_custom(SNR,d,M,modul_name,Ts)
    if nargin < 5
        Ts = 1;
    end
    data_rates = (1-d./(2.^d-1))'.*log2(M)./Ts;%bits/sec
    err = zeros(length(M),length(d),length(modul_name),3);
    for m = 1:length(M)
        fprintf('M = %d\n',M(m))
        for j = 1:length(d) 
            fprintf('d = %d\n',d(j))
            for i = 1:length(modul_name) 
                [err(m,j,i,1),err(m,j,i,2),err(m,j,i,3)] = hamming_sim(SNR,M(m),modul_name{i},d(j),54,3e7);%numofbits = 1e8
            end
        end
    end
