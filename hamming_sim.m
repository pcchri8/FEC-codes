function [biterr_c,blockerr_c,numofbits]=...
                          hamming_sim(SNR,M,modul_name,d,bnum,max_mem)
    if nargin<5
        max_mem = 16e8;%160MB
        bnum = 1;
    elseif nargin<6
        max_mem = 16e8;%160MB
    end
    n = 2^d-1;
    k = n-d;
    %% Memory management
    % calculating packet length so that the script (not the whole app) 
    % doesnt utilize more than max_mem bytes of RAM 
    % we need 2 copies of the data and some calculations
    % add a copy temporarily (fail safe side estimation)
    min_packet_sz = lcm(k,log2(M));%for the coding scheme to be feasible
    if min_packet_sz <= max_mem/16 && d*n <= max_mem/16
        cuts = 1;
        l = floor(max_mem/(4*min_packet_sz*4));
        data_sz = l*min_packet_sz;
    else
        error('Not enough memory.')
    end
    %% Init
    rng('default')
    H = hammgen(d);
    p = pb_err(M,modul_name,SNR);
    blockerr_c = 0;%block errors counter
    biterr_c = 0;%bit errors counter
    %Calculating syndrome matrix
    s = H';%Hamming only
    P =s(d+1:end,:);%P' so that data*P becomes (P'*data')' since data will be stored columnwise
    data = randi([0,1],data_sz,1,'single');%Generating data
    data = reshape(data,[],k);
%% Main part
    if d < 11
        parfor j = 1:bnum%if we want to use less memory but still a lot data
            temp = [mod(data*P,2),data];%Encoding
            bsc_ind = randperm(numel(temp),ceil(p*numel(temp)));%Simulating binary symmetric channel for large amount of data bit errors -> p*bits
            temp(bsc_ind) = mod(temp(bsc_ind)+1,2);%p*numel(bits) are swapped
            syn = mod(temp*s,2);%Calculating syndrome of noisy data
        %         temp = mod(temp + permute(sum(syn == permute(s,[3,2,1]),2) == d,[1,3,2]),2); %Faster for larger codewords  d>11
            for i =d+1:n%2^(n-k)-1
                ind = sum(syn == s(i,:),2);%ind is a size(data) array with values ranging [0,n-k] storing how many bits are equal in the comparison
                temp(ind==d,i) = mod(temp(ind==d,i) + 1,2);%hamming only(if n-k bits are equal the syndromes are the same so change the ith bit)
            end
            temp = temp(:,d+1:end);%Decoding (systematic code)
            biterr_c = biterr_c + sum(sum(~(temp==data)));%all the bits that are different are counted
            blockerr_c = blockerr_c + sum(logical(sum(~(temp==data),2)));%all the blocks that are different are counted
        end
    else
        parfor j = 1:bnum%if we want to use less memory but still a lot data
            temp = [mod(data*P,2),data];%Encoding
            bsc_ind = randperm(numel(temp),ceil(p*numel(temp)));%Simulating binary symmetric channel for large amount of data bit errors -> p*bits
            temp(bsc_ind) = mod(temp(bsc_ind)+1,2);%p*numel(bits) are swapped
            syn = mod(temp*s,2);%Calculating syndrome of noisy data
            temp = mod(temp + permute(sum(syn == permute(s,[3,2,1]),2) == d,[1,3,2]),2); %Faster for larger codewords  d>11
            temp = temp(:,d+1:end);%Decoding (systematic code)
            biterr_c = biterr_c + sum(sum(~(temp==data)));%all the bits that are different are counted
            blockerr_c = blockerr_c + sum(logical(sum(~(temp==data),2)));%all the blocks that are different are counted
        end
    end
    numofbits = data_sz*bnum*cuts;
    biterr_c = biterr_c/numofbits;
    blockerr_c = k*blockerr_c/numofbits;