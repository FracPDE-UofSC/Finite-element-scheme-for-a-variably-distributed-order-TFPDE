function result = function_g(a0,x,t)


% alpha = 0.2;
% result = t^2*sin(2*pi*x)
left= 0;  
right =1 ;
if  x == left
   result = t^(2-a0)*sin(2*pi*x);
elseif x==right
    result = t^(2-a0)*sin(2*pi*x);
else
   warining='error' 
end