function result= local_FE_basis_function_1D(x,vertices,basis_type,basis_index,der)

h = vertices(2)-vertices(1);
x_hat = (x-vertices(1))/h;
Jacobi = 1/h;

if basis_type == 101
    
    if der == 0
    result = reference_FE_basis_function_1D(x_hat,basis_type,basis_index,der);
    elseif der==1 
        result= Jacobi*reference_FE_basis_function_1D(x_hat,basis_type,basis_index,der);
    elseif der>1
        result = 0;
    else 
        warning='wrong of derivative of the local bais type 101';
        
    end
        
        
elseif basis_type==102
    if der==0
         result = reference_FE_basis_function_1D(x_hat,basis_type,basis_index,der);
    elseif der==1
         result = Jacobi*reference_FE_basis_function_1D(x_hat,basis_type,basis_index,der);
    elseif der==2
         result = Jacobi^2* reference_FE_basis_function_1D(x_hat,basis_type,basis_index,der);
    else
        warning = 'wrong of derivative of the local bais type 201';
    end
    
    
else
    warning='wrong of the basis_type';
end
