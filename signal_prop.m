%Task 1
%This script generates binary data, encodes them, simulates a binary
%symetric channel and then decodes them and calculates block error rate and
%ber.(Hamming only)
%
%Team #8
%Christof Petros 9928
%Grigoriou Stergios 9564
%Zacharioudaki Danai 9418
function [biterr_c,blockerr_c,numofbits]=...
                          signal_prop(SNR,M,modul_name,d,bnum,max_mem)
%     SNR = 20;
%     M = 1024;
%     modul_name = 'PAM';
%     bnum = 1;
%     max_mem = 6e9;
%     d = 4;
    if nargin<5
        max_mem = 4e9;%4GB
        bnum = 1;
    elseif nargin<6
        max_mem = 4e9;%4GB
    end
    n = 2^d-1;
    k = n-d;
    %% Memory management
    % calculating packet length so that the script (not the whole app) 
    % doesnt utilize more than max_mem bytes of RAM 
    % we need 2 copies of the data and some calculations
    % add a copy temporarily (fail safe side estimation)
    % + the generator matrix must be below 0.4 of the threshold as well
    min_packet_sz = lcm(k,M);%for the coding scheme to be feasible
    if min_packet_sz <= max_mem/32 && k^2 <=max_mem/80
        cuts = 1;
%         if max_mem > 4e9    
%             l = floor(4e9/(4*min_packet_sz*8));%If it is not needed we don't use more than 4GB
%             if l < 1
%                 l=1;
%             end
%         else
        l = floor(max_mem/(4*min_packet_sz*8));
%         end
        data_sz = l*min_packet_sz;
    elseif k^2 <= max_mem/64%Not too large generator matrix
        cuts = M;
        data_sz = k;
        l=1;
    else
        error('Not enough memory.')
    end
    %% Init
    rng('default')
    [H,G] = hammgen(d);
    p = pb_err(M,modul_name,SNR);
    blockerr_c = 0;%block errors counter
    biterr_c = 0;%bit errors counter
    %Calculating syndrome matrix
    e = sparse(eye(n));%Hamming only
    s = e*H';
%% Main part
    for cut = 1:cuts%if the block size is bigger than max_mem we cut it in parts
        for j = 1:bnum%if we want to use less memory but still a lot data
            fprintf('Batch %d out of %d\n',j,bnum)
            fprintf('Packets %d to %d\nPacket size %d blocks:\n',(j-1)*l+1,j*l,min_packet_sz/k)
            fprintf('Cut %d out of %d\n\n',cut,cuts)
            fprintf('Data Generation...')
            batchtime = tic;
            t = tic;
            data = randi([0,1],data_sz,1,'double');%Generating data
            data = reshape(data,k,[]);
            t = toc(t);
            fprintf('\t %.3f seconds\nEncoding...',t)
            t = tic;
            temp = mod(data'*G,2)';%Encoding matmul mod 2 
            t = toc(t);
            fprintf('\t\t\t %.3f seconds\nTransmitting...',t)
            t = tic;
            temp= bsc(temp,p);%Simulating channel
            t = toc(t);
            fprintf('\t\t %.3f seconds\nError correcting...',t)
            t = tic;
            syn = mod((H*temp)',2);%Calculating syndrome of noisy data
            for i =1:n%2^(n-k)-1
                ind = sum(syn == s(i,:),2);%ind is a size(data) array with values ranging [0,n-k] storing how many bits are equal in the comparison
                temp(i,ind==d) = mod(temp(i,ind==d) + 1,2);%hamming only(if n-k bits are equal the syndromes are the same so change the ith bit)
            end
            t = toc(t);
            fprintf('\t %.3f seconds\nDecoding...',t)
            t = tic;
            temp = temp(n-k+1:end,:);%Decoding (systematic code)
            t = toc(t);
            fprintf('\t\t\t %.3f seconds\nError calculations...',t)
            t = tic;
            %bter = sum(~(temp==data));
            biterr_c = biterr_c + sum(sum(~(temp==data)));%all the bits that are different are counted
            blockerr_c = blockerr_c + sum(logical(sum(~(temp==data))));%all the blocks that are different are counted
            t = toc(t);
            fprintf('%.3f seconds\n\nBatch total time:',t)
            t = toc(batchtime);
            fprintf('\t %.3f seconds\n\n\n',t)
        end
    end
    numofbits = data_sz*bnum*cuts;
    biterr_c = biterr_c/numofbits;
    blockerr_c = k*blockerr_c/numofbits;