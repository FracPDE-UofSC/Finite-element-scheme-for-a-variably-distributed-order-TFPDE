function [P,T]=generate_PT(left, right, h, basis_type)
num_of_elements = (right-left)/h; %num_of_elements
Nb = num_of_elements+1;%num of nodes
for i = 1:Nb
   x(i) = left+(i-1)*h; % info_of_mesh_nodes
end
if basis_type==101 % 1D linear fem 
  P = x;%每个节点对应的具体的x值
  T = zeros(2,num_of_elements);%每个单元对应的总体节点编号
  for i = 1:num_of_elements
      T(1,i) = i;
      T(2,i) = i+1;
  end
  
    
else %basis_type==102 % 1D quadratic fem  
    Nb = 2*num_of_elements+1;
   for i =1:Nb
    y(i) = left+(i-1)*h/2;
   end
   P = y;
   T = zeros(3,num_of_elements);
   for i = 1:num_of_elements
       T(1,i) = 2*i-1;
       T(2,i) = 2*i+1;
       T(3,i) = 2*i;
   end
   
end

end


    