%Task 1
%
%Team #8
%Christof Petros 9928
%Grigoriou Stergios 9564
%Zacharioudaki Danai 9418
close all
clear
clc
%% Init
rng('default')
SNR = 1;
M = 2;
modul_name = 'PAM';
% p = pb_err(M,modul_name,SNR);
p = 1;
[H,G,n,k] = hammgen(3);
% k = 4;
% n = 7;
l = 1e5;
data = randi(2,n*M*k*l,1) - 1;
data = reshape(data,k,[]);
cd_data = zeros(n,M*n*l);
for i = 1:n*M*l
    cd_data(:,i) = mod(data(:,i)'*G,2);
end
channel_data = bsc(cd_data,p);
% syn = mod((H*channel_data)',2);
% syn = binaryVectorToDecimal(syn);
% [e ,~] =  typical_lay(n,k,H);
% for i =1:2^(n-k)-1
%     channel_data(:,syn==i) = mod((channel_data(:,syn==i) + e(i+1,:)'),2);
% end
block_er_r = mod((H*channel_data),2); 
block_er_r = sum(logical(sum(block_er_r)))/(M*l*n);




