function result = func_bkn(t,k,n,tau,alpha)
 % t: a series of t: t1;t2;t3;;;;
 % k: index of the tk
 % n: current t at n
 % alpha: order of the fractional
 % tau: step of time
 % compute the binomial coefficient of b(k,n)

result = ( t(n)-t(k-1) )^(1-alpha)  - ( t(n)-t(k) )^(1-alpha);
result = result / ( t(k)-t(k-1) );
% result = result / tau;

