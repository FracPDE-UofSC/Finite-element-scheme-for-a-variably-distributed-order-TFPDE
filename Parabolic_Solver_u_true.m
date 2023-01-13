function [u]  = Parabolic_Solver_u_true(num_of_time)


% solve time fractional equation by using fem 1D case
% left right : domain
% h  : space size
% basis_type_trial
% basis_type_test

left  =  0; right = 1; %domain
Tstart = 0; Tfinal = 1; %time
boundary_type = -1; % Dirichlet boundary condition
%% check(alpha)
a0 = 0.1;  a1 = 0.5;
%% hat(alpha)
a00 = 0.6;  a11 = 0.8;
%%
r = 1;%Uniform Mesh
% r = 4/(a00);%Graded Mesh
num_of_elements = 32;
h = (right-left)/num_of_elements; % space st
Mm = num_of_time; % number of time step
dt = (Tfinal-Tstart)/(Mm); % time step
if r == 1
    disp('Uniform Mesh')
    for i = 1:Mm+1
        t(i) = ((i-1)/Mm )^r *(Tfinal-Tstart);
    end
else
    %% Graded mesh
%     r = 1/(a00);
    for i = 1:Mm+1
        t(i) = ((i-1)/Mm )^r *(Tfinal-Tstart);
    end
    disp('Graded Mesh')
end
% % t = sort(t);
basis_type = 101; % 1D linear fem
basis_type_trial=basis_type; % 1D
[P,T]=generate_PT(left, right, h, basis_type); % 1D linear fem
if basis_type_trial==101
    Pb_trial=P;
    Tb_trial=T;
else
    [Pb_trial,Tb_trial]=generate_PbTb(left, right, h, basis_type);
end
Pb = Pb_trial; Tb = Tb_trial;
% boundarynodes=generate_boundary_node_info(nbn,Nb,boundary_type)
Nb = size(Pb_trial,2);
boundarynodes=generate_boundary_node_info(left,right,Pb_trial,boundary_type);
matrix_size=compute_number_of_unknown(Nb);
% basis_type==101  101: 1D linear FE
% number_of_local_basis_trial=identify_number_of_local_basis(basis_type_trial);
number_of_local_basis_trial = size(Tb_trial,1) ;
number_of_local_basis_test =  size(Tb_trial,1) ;
A=assemble_matrix_1D('coe_fun_a',matrix_size,Tb_trial,Pb_trial,...
    number_of_local_basis_trial,basis_type,1,...
    number_of_local_basis_test,basis_type,1); % stiffness matrix
M = assemble_matrix_1D('coe_fun_m',matrix_size,Tb_trial,Pb_trial,...
    number_of_local_basis_trial,basis_type,0,...
    number_of_local_basis_test,basis_type,0); % Mass matrix

A_fixed = A;
M_fixed = M;
% initial solution

x0 = zeros(size(A_fixed,1),1); % initial vector x0
for i = 1:length(x0)
    x0(i) = func_u0(Pb(1,i)); % initial condition
end
x0_fixed = x0;
Uh = zeros(length(x0),length(t)); % store every time solution
Uh(:,1) = x0_fixed; %initial solution at start time t==0


%% Iterate time
for iTime = 1:length(t)-1
    oldLevel = iTime;    newLevel =  iTime+1;
    t1 = t(iTime+1);       t0 = t(iTime);   tau = t1 - t0;
    
    alpha1 = func_alpha(a0,a1,t1); %\check {alpha(t)}
    alpha2 = func_alpha(a00,a11,t1);%\hat {alpha(t)}
    interval = [alpha1,alpha2];% integration interval
    type = 106; %  6 Gauss_nodes
    type = 1020;% 20 Gauss_nodes
    [Gauss_weights,Gauss_nodes] = generate_Gauss(interval,type);
    val_int_w = func_int_w(a00,Gauss_weights,Gauss_nodes);
%     num_of_alpha =  ceil((alpha2-alpha1)/h );% number of alpha for  (\check(alpha()t),\hat(alpha(t)))
%         num_of_alpha =  ceil((alpha2-alpha1)/dt^(0.5));% number of alpha for  (\check(alpha()t),\hat(alpha(t)))
%     num_of_alpha = num_of_elements;
%     vector_ksi = linspace(alpha1,alpha2,num_of_alpha);
%     for j = 1:length(vector_ksi)-1
%         vector_alpha(j) = 1/2*(vector_ksi(j) + vector_ksi(j+1) );
%     end
    
    coef_bnn = func_coef_sum_bnn(a00,Gauss_weights,Gauss_nodes,t,tau,newLevel);
    coef_b2n = func_coef_sum_b2n(a00,Gauss_weights,Gauss_nodes,t,tau,newLevel);
    A_tilde =   coef_bnn* M_fixed/val_int_w + A_fixed ...
        + M_fixed/tau;
    b1 = assemble_vector_1D('func_f',a00,t1,alpha1,alpha2,matrix_size,Tb,Pb,number_of_local_basis_test, 101,0);
    
    sum_bm = func_vec_sum_bkn(x0,t,tau,a00,Gauss_weights,Gauss_nodes,Uh,newLevel);
    sum_b = M_fixed * sum_bm ;
    b_tilde = b1/val_int_w...
        - sum_b/val_int_w ...
        +  coef_b2n*M_fixed*Uh(:,1)/val_int_w...
        + M_fixed*Uh(:,oldLevel)/tau;
    %% treat bounadry condition
    [AA,b]=treat_Dirichlet_boundary(A_tilde,b_tilde,t1,'function_g',a00,Pb,boundarynodes);
    %% solve
%     x1=AA\b;
    x1 = cgs(AA,  b,1e-6 ,200);
    x0 = x1; % renew x0
    Uh(:,iTime+1) = x0;
end % of m for time steps




u = x0; % numerical solution at the terminal time
