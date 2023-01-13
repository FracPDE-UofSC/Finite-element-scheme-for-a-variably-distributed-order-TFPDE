function result =  func_sum_bm(func_sum_b,func_w,a0,func_bkn,vector_alpha,...
        t,newLevel,tau,Uh,x0)
% compute vector sum_bmn  k=2:n-1

[row,~] = size(Uh);
sum_bm = zeros(row,1);
rou = vector_alpha(2) - vector_alpha(1);
for i  = 1:length(vector_alpha)
    sum_bm  = sum_bm + rou/gamma(2-vector_alpha(i))...
        *feval ( func_w,a0,vector_alpha(i) )...
        *feval ( func_sum_b,x0,t,newLevel,vector_alpha(i),tau,Uh );
end
result = sum_bm;