function result = func_coef_sum_bnn(a00,Gauss_weights,Gauss_nodes,t,tau,n)


result = 0;
% rou = Gauss_nodes(2)-Gauss_nodes(1);
Nalpha =  length(Gauss_nodes);
for j = 1:Nalpha
   result = result  ...
       + Gauss_weights(j)...
        *func_w(a00,Gauss_nodes(j)) * 1/gamma(2-Gauss_nodes(j)) ...
       *func_bkn(t,n,n,tau,Gauss_nodes(j));
end
result = result;