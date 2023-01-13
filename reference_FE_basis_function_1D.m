function result= reference_FE_basis_function_1D(x_hat,basis_type,basis_index,der)



if basis_type == 101 %linear finite element in 1D
    if basis_index == 1
        if der == 0
            result = 1-x_hat;
        elseif der== 1
            result = -1;
        elseif der>1
            result = 0;
        else
            waring='wrong of the derivative of basis  of 101 of index 1';
        end
    elseif basis_index==2
        if der == 0
            result = x_hat;
        elseif der==1
            result = 1;
        elseif der>1
            result = 0;
        else
            waring='wrong of the derivative of basis  of 101 of index 2';
        end
    else
        waning = 'warong of basis_index of 101';
    end
    
    
    
elseif basis_type == 102%quadratic finite element in 1D
    if basis_index == 1
        if der==0
            result = 2*x_hat^2-3*x_hat+1;
        elseif der ==1
            result = 4*x_hat-3;
        elseif der==2
            result = 4;
        elseif der>2
            result =0;
        else
            waring = 'wrong of the derivative of basis  of 201 of index 1 ';
        end
    elseif basis_index==2
        if der==0
            result = 2*x_hat^2-x_hat;
        elseif der ==1
            result = 4*x_hat-1;
        elseif der==2
            result = 4;
        elseif der>2
            result =0;
        else
            waring = 'wrong of the derivative of basis  of 201 of index 2';
        end
    elseif basis_index==3
        if der==0
            result = -4*x_hat^2+4*x_hat;
        elseif der ==1
            result = -8*x_hat+4;
        elseif der==2
            result = -8;
        elseif der>2
            result =0;
        else
            waring = 'wrong of the derivative of basis  of 201 of index 3';
        end
    else
        warning = 'wrong of basis_index of 201';
    end
   
else
    warning='wrong of basis_type'
end