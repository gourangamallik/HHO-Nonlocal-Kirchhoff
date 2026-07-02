function [uTF,d,nIter] = Newton_solver(idofs,A_global,AS_global, b_global, Prev_uTF,Prev_d)
% ----------------------------------------------------------------------
%      Solve Kirchhoff Problem Using Newton's Method
% ----------------------------------------------------------------------
%
% Author: Gouranga Mallik
%
% ----------------------------------------------------------------------
N=length(b_global); dofs=[idofs, N+1];
uTF=sparse(N,1);beta=[Prev_uTF;Prev_d];
Iter_Err=1; nIter=0; tic;TOL=10^(-10); Newton_Err=[];relax_param=1;
while(Iter_Err>TOL && nIter<50)
    nIter=nIter+1;
    d = Prev_uTF'*(A_global*Prev_uTF); Md=1+d;
    Jc=2*A_global*Prev_uTF;
    F=[Md*AS_global*Prev_uTF-b_global; d-Prev_d];
    JA=Md*AS_global; Jb=AS_global*Prev_uTF; 
    Jacobian=[JA   Jb;
              Jc'  -1];

    beta(dofs) = beta(dofs) - relax_param * Jacobian(dofs,dofs)\F(dofs);
    uTF = beta(1:N); d = beta(N+1);

    Iter_Err=sqrt((uTF-Prev_uTF)'*(AS_global*(uTF-Prev_uTF)))/(uTF'*(AS_global*uTF));
    Prev_uTF=uTF; Prev_d=d;
    Newton_Err=[Newton_Err, Iter_Err];
end
%Newton_OC = compute_Newton_Convergence(Newton_Err,relax_param); %compute the rate of convergence of the Newton Iterations
end

