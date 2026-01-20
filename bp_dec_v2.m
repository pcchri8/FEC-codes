%Team #8
%Christof Petros 9928
%Grigoriou Stergios 9564
%Zacharioudaki Danai 9418

%% Belief Propagation Decoder
% This function uses the product-sum algorithm to decode an array of
% codewords given a parity check matrix H.
%
% -IN-
% enc(binary array with NaN values) is an array storing the data to be
% decoded. Each row is a codeword. 
%
% H(binary array) is the parity check matrix
%
% epsilon(non negative scalar) is the erasure probability
%
% -OUT-
%
%corrected(binary array) is the decoded codeword

function corrected = bp_dec_v2(enc,H,epsilon)
    n = size(H,2);
    s = H';
    maxiter = ceil((0.1+epsilon)*n)+10;%How many iterations of belief propagation
    corrected = zeros(size(enc));%To store the output
    tic
    parfor k = 1:size(enc,1)%For each codeblock in input data
        dec = enc(k,:);
        disp(k)
        temp = dec;
        for i = 1:maxiter
            erasures = isnan(dec);%Calculate the position of the erasures in codeblock
            ind = find(erasures);
            if isempty(ind)|| i == maxiter%If no erasures return word
                corrected(k,:) = dec;
                break
            elseif isequal(temp,dec) && i ~=1
                corrected(k,:) = dec;
                break
            else
                temp = dec;
            end
            %else
            for j = 1:length(ind)%for each erasure in the word
                other = isnan(dec);%Recalculate after each index
                other(ind(j)) = false;%Set current erasure out
                dec(ind(j)) = 0;%Destroy the edges of the erasure
                dec(other) = n+1;%Nan handle for the other erasures
                check_nodes = dec*s;%Calculate check_nodes output
                % (if for some node the output is >n it means that there 
                % was an erasure)
                dec(other) = nan;%Restore nan value 
                check_nodes(check_nodes>n) = nan;%Exploiting the nan handle
                check_nodes = mod(check_nodes,2);%Value now is 0 1 or NaN
                decision = check_nodes(~isnan(check_nodes)&s(ind(j),:));%
                % Keeping only the check nodes that are connected to our 
                % erasure and their ouput is not an erasure. 
                var_node = sum(decision);% They will be either all 1 or all 0
                if var_node > 0 
                    dec(ind(j)) = 1;
                elseif var_node == 0 && ~isempty(decision)
                    dec(ind(j)) = 0;
                else%if decision is empty
                    dec(ind(j)) = nan;
                end
            end
        end
    end
    toc