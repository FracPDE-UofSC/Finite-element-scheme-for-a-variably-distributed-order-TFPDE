function boundarynodes=generate_boundary_node_info(left,right,Pb_trial,boundary_type)
% nbn: number_of_boundary_nodes
% Nb : number_of_basis_functions
% boundary_type -1   1d Dirichlet
% boundary_type -2   1d Neuman
% boundary_type -3   1d Robin
domain(1) =left;
domain(2) = right;
nbn = length(domain);
boundarynodes = zeros(2,nbn);
if boundary_type == -1
    boundarynodes(1,:) = -1;
%     boundarynodes(2,1) = 1;
%     boundarynodes(2,2) = Nb;
    boundarynodes(2,:)=find(Pb_trial==left|Pb_trial==right);
elseif boundary_type == -2
    boundarynodes(1,:) = -2;
    boundarynodes(2,1) = 1;
    boundarynodes(2,2) = Nb;
else boundary_type == -3
    boundarynodes(1,:) = -3;
    boundarynodes(2,1) = 1;
    boundarynodes(2,2) = Nb;
end
end