function [dt,error_max]  = Parabolic_Solver_h_order(num_of_elements)


% solve time fractional equation by using fem 1D case
% left right : domain
% h  : space size
% basis_type_trial
% basis_type_test

left  =  0; right = 1; %domain
Tstart = 0; Tfinal = 1; %time
boundary_type = -1; % Dirichlet boundary condition
% num_of_elements = 2^4;
alpha1 = 0.5;
alpha2 = 1;
% num_of_element = ceil( num_of_elements^(0.5));
h = (right-left)/num_of_elements; % space st
Mm = num_of_elements; % number of time step
num_of_alpha = num_of_elements;
dt = (Tfinal-Tstart)/(Mm); % time step
t = linspace(Tstart,Tfinal,Mm); % discrete time
dt = t(2)-t(1);
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
A_fixed=assemble_matrix_1D('coe_fun_a',matrix_size,Tb_trial,Pb_trial,...
    number_of_local_basis_trial,basis_type,1,...
    number_of_local_basis_test,basis_type,1); % stiffness matrix
M_fixed = assemble_matrix_1D('coe_fun_m',matrix_size,Tb_trial,Pb_trial,...
    number_of_local_basis_trial,basis_type,0,...
    number_of_local_basis_test,basis_type,0); % Mass matrix
% initial solution
x0 = zeros(size(A_fixed,1),1); % initial vector x0
for i = 1:length(x0)
    x0(i) = func_u0(Pb(1,i),0); % initial condition
end
x0_fixed = x0;
Uh = zeros(length(x0),length(t)); % store every time solution
Uh(:,1) = x0_fixed; %initial solution at start time t==0

vector_ksi = linspace(alpha1,alpha2, num_of_alpha);
for j = 1:length(vector_ksi)-1
    vector_alpha(j) = 1/2*(vector_ksi(j) + vector_ksi(j+1) );
end
%% Iterate time
for iTime = 1:length(t)-1
    oldLevel = iTime;    newLevel =  iTime+1;
    t1 = t(iTime+1);       t0 = t(iTime);   tau = t1 - t0;

    coef_bnn = func_coef_sum_bnn(vector_alpha,t,tau,newLevel);
    coef_b2n = func_coef_sum_b2n(vector_alpha,t,tau,newLevel);
    A_tilde =   coef_bnn* M_fixed + A_fixed ;
    b1 = assemble_vector_1D('func_f',t1,alpha1,alpha2,matrix_size,Tb,Pb,number_of_local_basis_test, basis_type,0);
    
    sum_bm = func_vec_sum_bkn(x0,t,tau,vector_alpha,Uh,newLevel);
    sum_b = M_fixed * sum_bm ;
    b_tilde = b1...
        - sum_b ...
        +  coef_b2n*M_fixed*Uh(:,1);
    %% treat bounadry condition
    [AA,b]=treat_Dirichlet_boundary(A_tilde,b_tilde,t1,'function_g',Pb,boundarynodes);
    %% solve
    x1=AA\b;
    x0 = x1; % renew x0
    Uh(:,iTime+1) = x0;
end % of m for time steps



uh = x0; % numerical solution at the terminal time
u = zeros(size(uh));
for i = 1:length(u)
    u(i) = u_true(Pb(1,i),Tfinal); % true solution at the terminal time
end
res = u-uh;
error_max = max(abs(res));
kanU = [u,uh,res];
%% plot
% plot(Pb,u,'r') % plot true solution at the terminal time
% hold on
% plot(Pb,uh,'b--') % plot numerical solution at the terminal time
% legend ('u','uh')
