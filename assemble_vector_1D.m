function b = assemble_vector_1D(f,a0,t,alpha1,alpha2,matrix_size,Tb_trial,Pb_trial,number_of_local_basis_test,basis_type_test,der_test)

% display begin_of_assemble_vector_1D;
b = zeros(matrix_size(1),1);
Tb_test = Tb_trial;
Tb_trial = Tb_trial;
number_of_elements=size(Tb_trial,2);
for n=1:number_of_elements
    vertices=Pb_trial(Tb_trial(:,n));
    [Gauss_weights,Gauss_nodes]=generate_Gauss(vertices,106);
    for beta=1:number_of_local_basis_test
        int_value=Gauss_quadrature_volume_1D_f(f,a0,t,alpha1,alpha2, Gauss_weights, Gauss_nodes, vertices,...
            101,1,0,...
            basis_type_test,beta,der_test);
        b(Tb_test(beta,n))=b(Tb_test(beta,n)) + int_value;
    end
end

end