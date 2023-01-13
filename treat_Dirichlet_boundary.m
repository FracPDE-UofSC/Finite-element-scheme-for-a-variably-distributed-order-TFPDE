function [A,b]=treat_Dirichlet_boundary(A,b,t,function_g,a0,Pb,boundarynodes)
%-1: Dirichlet

% display begin_of_treat_Dirichlet_boundary;
[~,nbn] = size(boundarynodes);
for k=1:nbn
    
    if boundarynodes(1,k)==-1% type of the boundary condition
        
        i = boundarynodes(2,k);%
        
        A(i,:)=0;
        A(i,i)=1;
        b(i,1)=feval('function_g',a0, Pb(1,i),t);
        
        end
        
        
    end
end