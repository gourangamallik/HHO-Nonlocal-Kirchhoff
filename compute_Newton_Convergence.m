function [Newton_OC] = compute_Newton_Convergence(Newton_Err,relax_param)
% ----------------------------------------------------------------------
%      Compute the Convergence Rate for Newton's Iterations
% ----------------------------------------------------------------------
%
% Author: Gouranga Mallik
%
% ----------------------------------------------------------------------
J=length(Newton_Err);Newton_OC=zeros(1,J-2); %Rate of convergence of Newton iterations:
fprintf('Newton_iteration error and convergence rate for relaxation (default: 1) : %2.2f\n',relax_param)
for j=1:J
    if(j>=3)
        Newton_OC(j-2)=log(Newton_Err(j)/Newton_Err(j-1))/log(Newton_Err(j-1)/Newton_Err(j-2));
        Newton_OC(j-2)=log(Newton_Err(j)/Newton_Err(j-1))/log(Newton_Err(j-1)/Newton_Err(j-2));
        fprintf('%.3e   %2.2f\n', Newton_Err(j), Newton_OC(j-2));
    else
        fprintf('%.3e\n',Newton_Err(j));
    end
end 