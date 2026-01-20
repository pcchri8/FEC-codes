%this function approximates the values of xL and xR that need to be
%enquired for computation of the L and R parameters of an irregular LDPC
%code and thus optimally determine an approximation of the code rate for a
%given value of codeword length (n)
%lamda and ro (column vectors) are the optimal values computed through the optimization of
%the code rate that fits a system's requirements (we'll compute them in
%another function)

function [xL_approx, xR_approx]= approx_method_ldpc(n, lamda, ro)

    %get l_max and r_max
    l_max = size(lamda,1);
    r_max = size(ro,1);

    %initialize values
    L_approx = zeros(l_max,1);
    L = zeros(l_max,1);
    P_approx = zeros(r_max,1);
    P = zeros(r_max,1);
    
    %calculating L and L_approx
    %L(i) represends the number of parent nodes with i branches and
    %L_approx(i) is the approximation of L(i) (floor)
    i = 2:l_max;
    i = i';
    w_lamda = lamda(i)./i;
    L_approx(i) = floor(((w_lamda(i-1)./sum(w_lamda)).*n));
    L_approx(1) = floor(((lamda(1)./sum(w_lamda)).*n));
    L(i) = ((w_lamda(i-1)./sum(w_lamda)).*n);
    L(1) = ((lamda(1)./sum(w_lamda)).*n);
    
    %calculating P and P_approx
    %P(i) represends the number of children nodes with i branches and
    %P_approx(i) is the approximation of P(i) floor(ceil)
    j = 2:r_max;
    j = j';
    w_ro = ro(j)./j;
    P_approx(j) = floor(((w_ro(j-1)./sum(w_lamda)).*n));
    P_approx(1) = floor(((ro(1)./sum(w_lamda)).*n));
    P(j) = ((w_ro(j-1)./sum(w_lamda)).*n);
    P(1) = ((ro(1)./sum(w_lamda)).*n);
    
    %metric for approximation error which appears in the code rate
    %approximation error
    xR = P(j)-P_approx(j);
    A = sum(xR);
    
%     first approach didn't seem to work out, needed more custom approach 
%     f = [ones(l_max-1,1); zeros(r_max-1,1)];
%     intcon = l_max+r_max-2;
%     Aeq1 = [ones(1,l_max-1), zeros(1,r_max-1)];
%     Aeq2 = [i', -j'];
%     Aeq = [Aeq1; Aeq2];
%     beq1 = n - sum(L_approx(i));
%     beq2 = sum(j.*P_approx(j)) - sum(i.*L_approx(i));
%     beq = [beq1; beq2];

        %set the minimization problem
        xL_approx = optimvar('xl',l_max-1,'LowerBound',0,'UpperBound',1,'Type','integer');
        xR_approx = optimvar('xr',r_max-1,'LowerBound',0,'UpperBound',1,'Type','integer');
        prob = optimproblem('Objective',sum(xr));
        cons1 = sum(xl) == n - sum(L_approx(i));
        cons2 = sum(i.*xl) - sum(j.*xr) == sum(j.*P_approx(j)) - sum(i.*L_approx(i));
        cons3 = sum(xr) >= ceil(A);
        prob.Constraints.cons1 = cons1;
        prob.Constraints.cons2 = cons2;
        prob.Constraints.cons3 = cons3;
        problem = prob2struct(prob);
        %solve the minimization problem
        [sol,fval1,exitflag1,output1] = intlinprog(problem);
        
        %set the maximization problem (equivalent to minimization of the
        %opposite value)
        xL_approx = optimvar('xl',l_max-1,'LowerBound',0,'UpperBound',1,'Type','integer');
        xR_approx = optimvar('xr',r_max-1,'LowerBound',0,'UpperBound',1,'Type','integer');
        prob = optimproblem('Objective',-sum(xr));
        cons1 = sum(xl) == n - sum(L_approx(i));
        cons2 = sum(i.*xl) - sum(j.*xr) == sum(j.*P_approx(j)) - sum(i.*L_approx(i));
        cons3 = sum(xr) <= floor(A);
        prob.Constraints.cons1 = cons1;
        prob.Constraints.cons2 = cons2;
        prob.Constraints.cons3 = cons3;
        problem = prob2struct(prob);
        %solve the maximization problem
        [sol,fval2,exitflag2,output2] = intlinprog(problem);
        
        %from the two solutions we'll keep the one with MMSE 
        
end