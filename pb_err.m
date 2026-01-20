%Team #8
%Christof Petros 9928
%Grigoriou Stergios 9564
%Zacharioudaki Danai 9418

%% Bit error probability
% This function returns the bit error probability for a M-ary PAM,PSK or
% QAM modulation. Assuming gray coding and reasonable symbol SNR
% Pb = Ps/log2(M) | Pb bit error probability and Ps symbol error
% probability.

%-IN

% M is the order of the modulation (2^k)
% modul_name the name of the modulation accepted values are
% ('PAM','PSK','QAM')
% SNR is the average symbol SNR 

%-OUT

%p is the bit error probability of the modulation over AWGN 

function p = pb_err(M,modul_name,SNR)
    %SNRb = SNR/log2(M);%bit SNR
    k = log2(M);
    switch modul_name
        case 'PAM'
            p = (2*(M-1)/(M*k))*qfunc((6*SNR/(M^2-1))^.5);
        case 'PSK'
            p = 2*qfunc(sqrt(2*SNR)*sin(pi/M))/k;
        case 'QAM'
            p1 = 2*(M^.5-1)*qfunc((3*SNR/(M-1))^.5)/M^.5;
            p = (1-(1-p1)^2)/k;
        otherwise
            error('You entered unsupported modulation name!')
    end