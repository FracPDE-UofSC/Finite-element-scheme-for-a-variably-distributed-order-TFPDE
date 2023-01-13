function result = func_int_w(a00,Gauss_weights,Gauss_nodes)


result = 0;
% rou = Gauss_nodes(2)-Gauss_nodes(1);
Nalpha =  length(Gauss_nodes);
for j = 1:Nalpha
   result = result  ...
       + Gauss_weights(j)...
        *func_w(a00,Gauss_nodes(j)) ;
end
result = result;