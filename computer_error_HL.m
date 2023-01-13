function error = computer_error_HL(u_true,u_true_x,u,uh,der_u,Pb,Tb,basis_type)

%% Compute error of the finite element solution
% u_true:
% u_true_x:   diff(u_true,x)
% u: vector of the true solution
% uh: vector of the numerical solution
% der_u: derivative of the u_true   der_u==0  L2 norm; der_u==1 H1 semi-norm
% Pb: information of nodes
% Tb: infomation of elements
% ZhiweiYang 04/15/2019
%%
N = size(Tb,2); % Total elements = Tb 的列数
number_of_loccal_basis = size(Tb,1);% = Tb 的行数
error = 0;
I = 0;
for n = 1:N % elements from 1 to N
    vertices=Pb(Tb(:,n));% vertices of the ith_element
    [Gauss_weights,Gauss_nodes]=generate_Gauss(vertices,106);
    % Compute 每个单元的积分值 I
    Gpn=length(Gauss_nodes);
    I=0;
    for num_Gpn=1:Gpn
       %% 先算出每个积分点的函数值
        % 计算 us （真解）的在高斯积分点的函数值
        if der_u==0
            us = feval(u_true,Gauss_nodes(num_Gpn)); %计算真解的函数在高斯点的值
        elseif der_u == 1
            us = feval(u_true_x,Gauss_nodes(num_Gpn));%计算真解的导函数在高斯点的值
        else 
            warning='wrong input for derivative!'
        end
        % 计算 数值解在在高斯积分点的函数值
        value_of_uh = 0;
        for k = 1:number_of_loccal_basis
            der_test = der_u; % 计算几阶导数的函数值
            value_of_uh = value_of_uh +  ...
                             uh(Tb(k,n))*local_FE_basis_function_1D(Gauss_nodes(num_Gpn),vertices,basis_type,k,der_test);
        end
        value_of_int_function = (us - value_of_uh)^2;
        I = I + Gauss_weights(num_Gpn) * value_of_int_function;
    end
    
    error = error + I;
    I = 0;
end
error = sqrt(error);

