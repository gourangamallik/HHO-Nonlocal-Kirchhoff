function [A_G_global,A_GS_global,L2Mass_global,L2Mx1_global,b_global,Bd_dofs]= global_A_b(hho, A_G_local,A_GS_local,M_L2,Mx1_L2,b)
% ----------------------------------------------------------------------
%                   Create Global Matrices
% ----------------------------------------------------------------------
%
% Author: Gouranga Mallik
%
% ----------------------------------------------------------------------
ncell_dofs = hho.mesh.ncells * hho.elements{1}.ncell_dofs;
nedge_dofs = hho.mesh.nedges * hho.elements{1}.nedge_dofs;
ntotal_dofs = ncell_dofs + nedge_dofs;
A_G_global=sparse(ntotal_dofs,ntotal_dofs); A_GS_global=sparse(ntotal_dofs,ntotal_dofs); b_global=sparse(ntotal_dofs,1);
L2Mass_global=sparse(ncell_dofs,ncell_dofs); L2Mx1_global=zeros(1,ncell_dofs);
Bd_dofs = []; % To collect boundary edge dofs
for cell_i=1:hho.mesh.ncells
    element = hho.elements{cell_i};
    g_cell_index=(cell_i-1)*element.ncell_dofs+(1:element.ncell_dofs)';
    g_edge_index=[];
    for edge_i = 1:element.nedges
        i_global = hho.mesh.cell_edges{cell_i}(edge_i);
        if hho.mesh.cell_neighbors{cell_i}(edge_i)==0 % To determine boundary edges
            ipos = (i_global - 1) * element.nedge_dofs + (1:element.nedge_dofs)';
            Bd_dofs = [Bd_dofs; ipos + ncell_dofs];
        end
        ipos = (i_global - 1) * element.nedge_dofs + (1:element.nedge_dofs)';
        g_edge_index=[g_edge_index;ipos + ncell_dofs];      
    end
    L2Mass_global(g_cell_index,g_cell_index) = L2Mass_global(g_cell_index,g_cell_index) + M_L2{cell_i};
    L2Mx1_global(g_cell_index) = L2Mx1_global(g_cell_index) + Mx1_L2{cell_i};
    g_index=[g_cell_index;g_edge_index];
    A_G_global(g_index,g_index)=A_G_global(g_index,g_index)+A_G_local{cell_i};
    A_GS_global(g_index,g_index)=A_GS_global(g_index,g_index)+A_GS_local{cell_i};
    b_global(g_cell_index)=b_global(g_cell_index)+b{cell_i};
end
end


