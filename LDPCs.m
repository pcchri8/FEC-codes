%Demo for how to use the suite
clear
clc
close
rng('default')
%use this to save time
data = load("li_ri_50_50.mat");%or load("li_ri_50_50.mat");
li = data.li_opt{1};
ri = data.ri_opt{1};
opt_rate = data.rate_opt;
clear data
%or use this [ri_opt,li_opt,rate_opt,rx,lx,rate] = li_ri_opt(50,50,100);
%li = li_opt{1};
%ri = ri_opt{1}
%opt_rate = rate_opt;
n = 1024;
epsilon = 1/7;
%The next code is how to find a good n
% newratearr = zeros(90,1);
%  for   n = 30:100
% %[Rf,Lf,newrate] = Li_Ri_approx_opt(ri,li,n);
%     [Rf,Lf,newratearr(n-29)] = Li_Ri_approx_opt_v1_2(ri,li,n);
%  end
%  n = 29+find(newratearr == max(newratearr));
[Rf,Lf,newrate] = Li_Ri_approx_opt_v1_2(ri,li,n);
rmax = find(Rf>0,1,'last');
lmax = find(Lf>0,1,'last');
Rf = Rf(1:rmax);
Lf = Lf(1:lmax);
n= sum(round(Lf));
H = poly2mat(Rf,Lf);
% H = hammgen(3);
% n = 7;
data = zeros(100,n);%encoded data
channel_data = data;%simulating bec
add_er = rand(100,n);
channel_data(add_er<epsilon) = nan;
corrected = bp_dec(channel_data,H,epsilon);
%Calculating Performance
iscor = zeros(size(data,1),1);
corr_bits = iscor;
for i = 1:size(data,1)
    corr_bits(i) = sum(corrected(i,:)==data(i,:));
    iscor(i) = corr_bits(i)==n;
end
block_err = 1-sum(iscor)/length(iscor);
ber = 1-sum(corr_bits)/numel(data);

