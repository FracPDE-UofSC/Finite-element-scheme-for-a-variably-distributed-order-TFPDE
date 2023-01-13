function result = func_w(a0,alpha)

% if 1-a0==0
%     result = 0;
% else
result = gamma(2-a0-alpha+1)/( gamma(3-a0)  );
% result = 1;
end