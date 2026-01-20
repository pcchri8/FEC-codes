%Team #8
%Christof Petros 9928
%Grigoriou Stergios 9564
%Zacharioudaki Danai 9418

%% Polynomials to Parity Check Matrix (Tanner Graph)
% This function given the integer Lambda and Rho coefficients calculates
% the parity check matrix for the corresponding code.
%
% -IN-
%
% Rf,Lf(vectors) the coefficients of the Rho and Lambda polynomials in ascending
% degree order
%
%-OUT-
%
% H is the parity check matrix
function H = poly2mat(Rf,Lf)
    %% Initializing
    L = round(Lf);%making sure they are integers
    R = round(Rf);
    n = sum(L);
    k = n - sum(R);
    H = zeros(n-k,n);
    row_deg = [];
    column_deg = [];
    %% Storing the degree nodes
    %Variable nodes
    for i = 1:length(L)
        if L(i) ~= 0
            row_deg = [row_deg;repmat(i,[L(i),1])];
        end
    end
    %Check nodes
    for i = 1:length(R)
        if R(i) ~= 0
            column_deg = [column_deg;repmat(i,[R(i),1])];
        end
    end
    %Suffling them (uniform distribution)
    row_deg = row_deg(randperm(n));
    column_deg = column_deg(randperm(n-k));
    %Suffling the edges 
    sigma = randperm(sum(L.*(1:length(L))'));
    j = 1;
    for i = 1:length(sigma)
        if row_deg(j)>0%If the variable node has remaining edges
            col_ind = sum(cumsum(column_deg)<sigma(i))+1;%Calculate the corresponding check node
            H(col_ind,j) = H(col_ind,j) + 1;%Connect it to the corresponding check node
            row_deg(j) = row_deg(j) - 1;%Update variable node's remaining edges
        else 
            j = j + 1;%Go to the next variable node
            col_ind = sum(cumsum(column_deg)<sigma(i))+1;
            H(col_ind,j) = H(col_ind,j) + 1;
            row_deg(j) = row_deg(j) - 1;
        end
    end
    H = mod(H,2);