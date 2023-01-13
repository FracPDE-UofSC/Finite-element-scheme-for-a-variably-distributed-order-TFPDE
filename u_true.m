function result = u_true(a0,x,t)

% result = exp(x)*sin(x);
% result = exp(x);


result = t^(2-a0)*sin(2*pi*x);
% result = t^4*sin(2*pi*x);
% 
% result  =  x^4*(x-1)*t^(3/2);
% result  = ( t^(2+alpha) +t +2 )* sin(pi*x);