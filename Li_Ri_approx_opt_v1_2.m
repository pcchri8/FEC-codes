%Team #8
%Christof Petros 9928
%Grigoriou Stergios 9564
%Zacharioudaki Danai 9418

%% Estimating the optimal integer Lambda_i and Rho_i for a given codeword length n
% This function stimates the optimal integer Lambda_i and Rho_i given the
% real lambda_i and rho_i and a codeword length n. It then returns them
% along with the new design_rate of the code.
%
% -IN-
% ri,li(vectors) like output of li_ri_opt
%
% n(int scalar)
%
% -OUT-
% Rf(vector) the optimal Rho_i polynomial coefficients in 
% ascending degree order
%
% Lf(vector) the optimal Lambda_i polynomial coefficients in 
% ascending degree order
%
% newrate(scalar) the design rate of the optimal Lf and Rf coefficients
%
%
function [Rf,Lf,newrate] = Li_Ri_approx_opt_v1_2(ri,li,n)
    lmax = length(li);
    rmax = length(ri);
    Li = li(2:end)./(2:lmax)';
    Ri = ri(2:end)./(2:rmax)';
    Ri = (Ri/sum(Li))*n;
    Li = (Li/sum(Li))*n;
    xr = Ri-floor(Ri);
    A = sum(xr);
    opts=optimoptions(@intlinprog,'Display','off');
    prob1 = optimproblem('ObjectiveSense','minimize');
    xl_hat = optimvar('xl_hat',lmax-1,'Type','integer','LowerBound',0,'UpperBound',1);
    xr_hat = optimvar('xr_hat',rmax-1,'Type','integer','LowerBound',0,'UpperBound',1);
    prob1.Objective = sum(xr_hat);
    prob1.Constraints.cons1 = sum(xl_hat) == n - sum(floor(Li));
    prob1.Constraints.cons2 = sum(xl_hat.*(2:lmax)') - sum(xr_hat.*(2:rmax)') ==...
        sum(floor(Ri).*(2:rmax)') -sum(floor(Li).*(2:lmax)');
    prob1.Constraints.cons3 = sum(xr_hat) >= ceil(A);
    l1 = solve(prob1,'Options',opts);
    
    prob2 = optimproblem('ObjectiveSense','minimize');
    xl_hat = optimvar('xl_hat',lmax-1,'Type','integer','LowerBound',0,'UpperBound',1);
    xr_hat = optimvar('xr_hat',rmax-1,'Type','integer','LowerBound',0,'UpperBound',1);
    prob2.Objective = -sum(xr_hat);
    prob2.Constraints.cons1 = sum(xl_hat) == n - sum(floor(Li));
    prob2.Constraints.cons2 = sum(xl_hat.*(2:lmax)') - sum(xr_hat.*(2:rmax)') ==...
        sum(floor(Ri).*(2:rmax)') -sum(floor(Li).*(2:lmax)');
    prob2.Constraints.cons3 = sum(xr_hat) <= floor(A);
    l2 = solve(prob2,'Options',opts);
    
    if abs(sum(l1.xr_hat) - A) < abs(sum(l2.xr_hat) - A)
        l = l1;
    elseif abs(sum(l1.xr_hat) - A) > abs(sum(l2.xr_hat) - A)
        l = l2;
    else
        R1 = floor(Ri) + l1.xr_hat;
        R2 = floor(Ri) + l2.xr_hat;
        newrate1 = 1 - sum(R1)/n;
        newrate2 = 1 - sum(R2)/n;
        if newrate1>=newrate2
            l=l1;
        else
            l=l2;
        end
    end


    Rf = floor(Ri) + l.xr_hat;
    Lf = floor(Li) + l.xl_hat;
    newrate = 1 - sum(Rf)/n;
    Rf = [0;Rf];
    Lf = [0;Lf];