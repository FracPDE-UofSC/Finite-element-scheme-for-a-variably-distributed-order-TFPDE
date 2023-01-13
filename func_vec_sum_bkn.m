function result = func_vec_sum_bkn(x0,t,tau,a00,Gauss_weights,Gauss_nodes,Uh,newLevel)


% rou = Gauss_nodes(2)-(1);
Nalpha =  length(Gauss_nodes);
vec_sum_bkn = zeros(length(x0),1);
for j = 1:Nalpha
    vec = zeros(size(x0));
    for k = 2:newLevel-1
        vec = vec...
            + (func_bkn(t,k,newLevel,tau,Gauss_nodes(j)) - func_bkn(t,k+1,newLevel,tau,Gauss_nodes(j)))...
            *Uh(:,k);
    end
    vec_sum_bkn = vec_sum_bkn...
          + Gauss_weights(j)...
          *func_w(a00,Gauss_nodes(j)) * 1/gamma(2-Gauss_nodes(j))*vec;
end
result = vec_sum_bkn;