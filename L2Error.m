function abs_error = L2Error(hho,A_L2, u,u_exact,tag)
% ----------------------------------------------------------------------
%      Compute L2-Error
% ----------------------------------------------------------------------
%
% Author: Gouranga Mallik
%
% ----------------------------------------------------------------------
    u_interp = HHOInterpolate(hho, u_exact);
    u_diff = u - u_interp;
    
    abs_error=0;l2norm_u = 0;
    for cell_i = 1:hho.mesh.ncells
        element = hho.elements{cell_i};
        u_diff_local_val=zeros(element.ncell_dofs,1);u_local_val=zeros(element.ncell_dofs,1);
        u_diff_local_val(1:element.ncell_dofs)=u_diff((cell_i-1)*element.ncell_dofs+(1:element.ncell_dofs));
        u_local_val(1:element.ncell_dofs)=u_interp((cell_i-1)*element.ncell_dofs+(1:element.ncell_dofs));
        abs_error = abs_error+u_diff_local_val'*(A_L2{cell_i}*u_diff_local_val);
        l2norm_u = l2norm_u + u_local_val'*(A_L2{cell_i}*u_local_val);
    end
    if(tag=='relative')
        abs_error=sqrt(abs_error)/sqrt(l2norm_u);
    else
        abs_error=sqrt(abs_error);
    end
    
end

