function [e ,s] =  typical_lay(n,k,h)
%     s = decimalToBinaryVector(0:2^(n-k)-1);
    a = zeros(2^(n-k),1);
    a(1) = 0;
    for i = 0:n-1
        a(i+2) = 2^i;
    end
    h = binaryVectorToDecimal(h);
    b = zeros(2^(n-k),3);
    for i = 1:n-k
        b(:,i) = bitand(a,h(i));
        temp = mod(sum(decimalToBinaryVector(b(:,i)),2),2);
        b(:,i) = temp;
    end
    e = decimalToBinaryVector(a);
    s = b;
    decs = binaryVectorToDecimal(s);
    [~,idx] = sort(decs);
    s = s(idx,:);
    e = e(idx,:);