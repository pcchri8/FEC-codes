%this function models a binary erasure channel with input arguments the
%desired bit erasure probability and the transmitted binary data vector and as output
%a binary data vector containg erasures at the bits that failed to be
%transmitted

function out_data= bec(pb_er, in_data)
    
    rn = rand(1, length(in_data));
    out_data = in_data;
    
    for i = 1:length(in_data)
        if rn(i) < pb_er
            out_data(i) = NaN; 
        end
    end
end
