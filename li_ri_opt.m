%Team #8
%Christof Petros 9928
%Grigoriou Stergios 9564
%Zacharioudaki Danai 9418

%% Calculating the optimal non integer lambda_i and rho_i
% This function calculates the optimal real lambda_i and rho_i for each
% combination of lambda_max and rho_max up to lmax and rmax. It then
% returns the optimal real lambda_i and rho_i over all the combinations
% based on the design rate along with an array of the design rates of the
% optimal lambda_i and rho_i for each combinations and those extra lambda_i 
% and rho_i. 
%
% -IN-
% rmax(int>=2) the maximum rho_max
%
% lmax(int>=2) the maximum lambda_max
%
% howmanyxs(int>=3: default=100) is the resolution of the discrete x
% values of the optimization
%
% epsilon is the erasure probability
%
% verbose(default = 0) for displaying the progress of the function
%
% -OUT-
% ri_opt(vector) the overall optimal rho_i polynomial coefficients in 
% ascending degree order
%
% li_opt(vector) the overall optimal lambda_i polynomial coefficients in 
% ascending degree order
%
% rate_opt(scalar) the design rate of overall the optimal li and ri 
% coefficients
%
% rx,lx (cell) containing the optimal coefficients for each 2:rmax and
% 2:lmax combination 
%
% rate(array) contains the  design rate of the optimal li and ri 
% coefficients for each 2:rmax and 2:lmax combination
%

function [ri_opt,li_opt,rate_opt,rx,lx,rate] = li_ri_opt(rmax,lmax,howmanyxs,epsilon,verbose)
    %% Setting Defaults
    if nargin<4
        epsilon=1/3;
        verbose=0;
    elseif nargin < 5
        verbose = 0;
    end
    %% Initializing
    ravg = zeros(rmax-1,lmax-1);%starting from 2
    r = ravg;
    rs = 2:rmax;
    x = 1/howmanyxs:1/howmanyxs:1;%discrete x values
    lx = cell(rmax-1,lmax-1);
    rx = cell(rmax-1,lmax-1);
    %% calculating ravg and r
    for ri = 1:rmax-1
        for j = 1:lmax-1
            ravg(ri,j) = mean(2:rs(ri));
            r(ri,j) = floor(ravg(ri,j));
        end
    end
    %% Calculating optimal ri and li for each rmax and lmax
    rate = zeros(rmax-1,lmax-1);
    opts=optimoptions(@linprog,'Display','off');
    for ri = 1:rmax-1 
        t= tic;
        for li = 1:lmax-1
            if verbose
                fprintf('rmax = %d, lmax = %d\n',ri+1,li+1)
            end
            %% Optimization setup
            prob = optimproblem('ObjectiveSense','maximize');%The max prob
            lambda = optimvar('lambda',li,'LowerBound',0);%lmax-1 because we want i >=2
            prob.Objective = sum(lambda./(2:li+1)');
            prob.Constraints.cons1 = sum(lambda) == 1;
            cons2 = optimconstr(howmanyxs);%discrete x values
            cons3 = (r(ri,li)-1)*r(ri,li)*(r(ri,li)+1-ravg(ri,li))/ravg(ri,li) +...
                r(ri,li)*(ravg(ri,li)- r(ri,li)*(r(ri,li)+1-ravg(ri,li)))/ravg(ri,li);
            for j = 1:howmanyxs
                summa = 0;
                for k = 1:li
                    summa = summa + lambda(k)*...
                        (1-r(ri)*(r(ri)+1-ravg(ri))/ravg(ri)*(1-x(j))^(r(ri)-1) - ...
                        (ravg(ri) - r(ri)*(r(ri)+1-ravg(ri)))/ravg(ri)*(1-x(j))^(r(ri)))^k;
                end
                cons2(j) = epsilon*summa - x(j) <=0;
            end
            prob.Constraints.cons2 = cons2;
            prob.Constraints.cons3 = lambda(1)<=1/(epsilon*cons3);
            l = solve(prob,'Options',opts);
            if isempty(l.lambda)
                rate(ri,li) = 0;
                lx{ri,li} = [];
                rx{ri,li} = [];
            else
                rhos = zeros(ri+1,1);
                rhos(r(ri)-1) = r(ri)*(r(ri)+1-ravg(ri))/ravg(ri);
                rhos(r(ri)) = (ravg(ri)-r(ri)*(r(ri)+1-ravg(ri)))/ravg(ri);
                temp = 1-sum(rhos./(1:ri+1)')/sum(l.lambda./(2:li+1)');
                rate(ri,li) = temp;
                lx{ri,li} = [0;l.lambda];
                rx{ri,li} = rhos;
            end
        end
        if verbose
            toc(t)
        end
    end
    rate_opt = max(rate,[],'all');
    li_opt = lx(rate==rate_opt);
    ri_opt = rx(rate==rate_opt);