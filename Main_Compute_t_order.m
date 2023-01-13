
clc;
clear;


number_of_time = [8,16,24,32]';
nSamples = length(number_of_time);
% [u]  = Parabolic_Solver_u_true(2^10);
% save u_data u

% load()
[t,error_max,order] = deal( zeros(size(number_of_time)) );
load('u_data.mat')
for iSample = 1:nSamples
%     tic
%     currentTime = number_of_time(iSample)
%     currentSample = sprintf(' iSample = %i',iSample)
    [t(iSample),error_max(iSample)] = Parabolic_Solver_t_order(number_of_time(iSample),u);
%     toc
end

%% compute error_order
for i = 2:nSamples
    HH = t(i-1); hh = t(i);
    error_H = error_max(i-1); error_h = error_max(i);
    order(i) = log(error_H/error_h) / ( log(HH/hh) );
end
t_order = order;
table(number_of_time,error_max,t_order)