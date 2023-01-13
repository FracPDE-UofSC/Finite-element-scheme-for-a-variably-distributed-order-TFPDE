function result=FE_basis_fun_1D(x,vertices,basis_type,basis_index,der)
% return local_basis_fun_1D
% x : symbolic
% basis_type; linear101 or quaratic 102
% basis_index : 1 ÓÒ°ëÖ§   2×ó°ëÖ§
% der£ºorder of derivative
%101: 1D linear FE
%% 101: 1D linear FE

h = vertices(2)-vertices(1);

if basis_type==101%101: 1D linear FE
    h = vertices(2)-vertices(1);
    %% %ÓÒ°ëÖ§
    if basis_index==1%ÓÒ°ëÖ§
        if der==0
            result=(vertices(2)-x)/h;
        elseif der==1
            result=-1/h;
        elseif der>1
            result=0;
        else
            warning='wrong input for derivative!'
        end
        %% ×ó°ëÖ§
    elseif basis_index==2%×ó°ëÖ§
        if der==0
            result=(x-vertices(1))/h;
        elseif der==1
            result=1/h;
        elseif der>1
            result=0;
        else
            warning='wrong input for index!'
        end
        %% quadratic fe
    end
elseif basis_type==102 %quadratic fe
    h = vertices(2)-vertices(1);
    if basis_index==1%ÓÒ°ëÖ§
        if der==0
            result=(vertices(2)-x)/h;
        elseif der==1
            result=-1/h;
        elseif der>1
            result=0;
        else
            warning='wrong input for derivative!'
        end
    elseif basis_index==2%×ó°ëÖ§
        if der==0
            result=(x-vertices(1))/h;
        elseif der==1
            result=1/h;
        elseif der>1
            result=0;
        else
            warning='wrong input for index!'
        end
    end
end
end




