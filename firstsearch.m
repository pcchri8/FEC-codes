clear
clc
close all
lmax = 16;
rmax = 32;
prob = optimproblem('ObjectiveSense','maximize');
lambda = optimvar('lambda',lmax-1,'LowerBound',0);
prob.Objective = sum(lambda./(2:lmax)');
prob.Constraints.cons1 = sum(lambda) == 1;
cons2 = optimconstr(100);
ravg = zeros(rmax-lmax,1);
r = ravg;
epsilon = 0.001;
rs = lmax+1:rmax;
x = 0.01:0.01:1;
for i = 1:rmax-lmax
    ravg(i) = mean(2:rs(i));
    r(i) = floor(ravg(i));
end
for i = 14%:rmax-lmax
    cons4 = (r(i)-1)*r(i)*(r(i)+1-ravg(i))/ravg(i) + r(i)*(ravg(i)...
        - r(i)*(r(i)+1-ravg(i)))/ravg(i);
    for j = 1:100
        summa = 0;
        for k = 1:lmax-1
            summa = summa + lambda(k)*...
                (1-r(i)*(r(i)+1-ravg(i))/ravg(i)*(1-x(j))^(r(i)-1) - ...
                (ravg(i) - r(i)*(r(i)+1-ravg(i)))/ravg(i)*(1-x(j))^(r(i)))^k;
        end
        cons2(j) = epsilon*summa - x(j) <=0;
    end
    
    prob.Constraints.cons2 = cons2;
    prob.Constraints.cons3 = lambda>=0;
    prob.Constraints.cons4 = lambda(1)<=1/(epsilon*cons4);
end
