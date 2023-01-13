function A=assemble_matrix_1D(coe_fun,matrix_size,Tb_trial,Pb_trial,...
    number_of_local_basis_trial,basis_type_trial,der_trial,...
    number_of_local_basis_test,basis_type_test,der_test)
%%
% display begin_of_assemble_matrix_1D;
A=sparse(matrix_size(1),matrix_size(2));
Tb_test = Tb_trial;
Tb_trial = Tb_trial;
number_of_elements=size(Tb_trial,2);
for n=1:number_of_elements
    vertices=Pb_trial(Tb_trial(:,n));
    [Gauss_weights,Gauss_nodes]=generate_Gauss(vertices,106);
    for alpha=1:number_of_local_basis_trial
        for beta=1:number_of_local_basis_test
            int_value=Gauss_quadrature_volume_1D(coe_fun, Gauss_weights, Gauss_nodes, vertices,...
                basis_type_trial,alpha,der_trial,...
                basis_type_test,beta,der_test);
            A(Tb_test(beta,n),Tb_trial(alpha,n))=A(Tb_test(beta,n),Tb_trial(alpha,n))+int_value;
        end
    end
    
end