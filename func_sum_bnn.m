function result = func_sum_bnn(func_w,func_bkn,vector_alpha,t,k,n,tau)
% compute coef of the matrix  coef*M*un
sum_bnn = 0;
rou = vector_alpha(2) - vector_alpha(1);
for i  = 1:length(vector_alpha)
    sum_bnn  = sum_bnn + rou/gamma(2-vector_alpha(i))...
        *feval ( func_w,vector_alpha(i) )...
        *feval ( func_bkn,t,k,n,tau,vector_alpha(i) );
end
result = sum_bnn;