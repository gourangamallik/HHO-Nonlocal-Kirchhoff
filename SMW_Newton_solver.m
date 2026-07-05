function [uTF,d,nIter] = SMW_Newton_solver(hho, AS_local, A_global,AS_global, b_global, Prev_uTF,Prev_d)
% ----------------------------------------------------------------------
%   Solve Kirchhoff Problem Using Sherman-Morrison-Woodbury Method
% ----------------------------------------------------------------------
%
% Author: Gouranga Mallik
%
% ----------------------------------------------------------------------
Iter_Err=1; nIter=0; tic;TOL=10^(-10); Newton_Err=[]; relax_param=1;
fprintf("Convergence of the Newton iterations:\n");
while(Iter_Err>TOL && nIter<50)
    nIter=nIter+1;
    d = Prev_uTF'*(A_global*Prev_uTF); Md=1+d;
    Jc=2*A_global*Prev_uTF;
    F=Md*AS_global*Prev_uTF-b_global;
    g=d-Prev_d;
    Jb=AS_global*Prev_uTF; 
    
    Md_A = cellfun(@(A) Md*A, AS_local, "UniformOutput",false); %Multiply each cell matrix by Md
    InvA_F=DiffSolver(hho, Md_A, F);
    InvA_b=DiffSolver(hho, Md_A, Jb);
    
    L=F+(g-Jc'*InvA_F)*Jb/(1+Jc'*InvA_b);
    uTF = Prev_uTF -  DiffSolver(hho, Md_A, L);
    d = Prev_d + (g-Jc'*InvA_F)/(1+Jc'*InvA_b);
    Iter_Err=sqrt((uTF-Prev_uTF)'*(AS_global*(uTF-Prev_uTF)))/(uTF'*(AS_global*uTF));
    Prev_uTF=uTF; Prev_d=d;
    Newton_Err=[Newton_Err, Iter_Err];
    fprintf('%d  %.4e \n',nIter,Iter_Err);
end
end

