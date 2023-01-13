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
N = size(Tb,2); % Total elements = Tb ������
number_of_loccal_basis = size(Tb,1);% = Tb ������
error = 0;
I = 0;
for n = 1:N % elements from 1 to N
    vertices=Pb(Tb(:,n));% vertices of the ith_element
    [Gauss_weights,Gauss_nodes]=generate_Gauss(vertices,106);
    % Compute ÿ����Ԫ�Ļ���ֵ I
    Gpn=length(Gauss_nodes);
    I=0;
    for num_Gpn=1:Gpn
       %% �����ÿ�����ֵ�ĺ���ֵ
        % ���� us ����⣩���ڸ�˹���ֵ�ĺ���ֵ
        if der_u==0
            us = feval(u_true,Gauss_nodes(num_Gpn)); %�������ĺ����ڸ�˹���ֵ
        elseif der_u == 1
            us = feval(u_true_x,Gauss_nodes(num_Gpn));%�������ĵ������ڸ�˹���ֵ
        else 
            warning='wrong input for derivative!'
        end
        % ���� ��ֵ�����ڸ�˹���ֵ�ĺ���ֵ
        value_of_uh = 0;
        for k = 1:number_of_loccal_basis
            der_test = der_u; % ���㼸�׵����ĺ���ֵ
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

