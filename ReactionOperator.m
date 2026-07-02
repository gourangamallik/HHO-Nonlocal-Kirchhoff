% ----------------------------------------------------------------------
%                        HHO Reaction Operator
% ----------------------------------------------------------------------
%
% Author: Gouranga Mallik
%
% ----------------------------------------------------------------------

function [MR_T,MRx1_T] = ReactionOperator(hho, cell_i)
%REACTION OPERATOR Compute the local HHO reaction term (u,v)_T
%   hho         The HHO data structure
%   cell_i      The index of the cell T
% Returns the local reaction operator matrix MRT

    element = hho.elements{cell_i};

    % Initialise matrices
    MR_T=zeros(element.ncell_dofs); 
    MRx1_T=zeros(1,element.ncell_dofs);
    % Compute quadrature over each edge of the cell
    for edge_i = 1:element.nedges
        % Compute volumetric terms (u,v)_T
        for iqn = 1:element.ncell_qnodes
            % The quadrature point (xQN,yQN) and the weight wQN
            xQN = element.cell_quad_rules{edge_i}(iqn,1);
            yQN = element.cell_quad_rules{edge_i}(iqn,2);
            wQN = element.cell_quad_rules{edge_i}(iqn,3);
                       % Compute cell value term
            for i = 1:element.ncell_dofs
               phi_i = element.cell_values_in_cell(edge_i, i, iqn);
               MRx1_T(i) = MRx1_T(i) + wQN * phi_i; % For (u,1)_T
               for j = 1:element.ncell_dofs
                  phi_j = element.cell_values_in_cell(edge_i, j, iqn);
                  MR_T(i, j) = MR_T(i, j) + wQN * phi_i * phi_j;
               end
            end
        end
    end
   
end

