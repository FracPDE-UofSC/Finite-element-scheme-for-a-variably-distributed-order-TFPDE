function result = func_sum_b(x0,t,newLevel,alpha,tau,Uh)


sum_b = zeros(size(x0));
for k = 2 : newLevel-1
    coef_bk = func_bkn(t,k,newLevel,tau,alpha)...
        - func_bkn(t,k+1,newLevel,tau,alpha);
    diff_u = Uh(:,k);
    sum_b = sum_b + coef_bk*diff_u ;
end % of k

result = sum_b;