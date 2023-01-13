function result=FE_function_1D(x,vertices,uh_local,basis_type,der)

number_of_local_basis=length(uh_local);

result=0;

for k=1:number_of_local_basis
    
    result=result+uh_local(k)*FE_basis_fun_1D(x,vertices,basis_type,k,der);
    
end

