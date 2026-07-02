% ----------------------------------------------------------------------
%   Main File to Run the Nonlocal Kirchhoff Problem
% ----------------------------------------------------------------------
%
% Author: Gouranga Mallik
%
% Last Modification: June-2026
% ----------------------------------------------------------------------

%% Nonlocal Kirchoff Problem: 
%
%         { -(1 + \|\nabla u\|^2 )\Delta u = f     in Omega
%         { u = 0     on Boundary(Omega)
%
close all;format long
% -------------------------------------------------------------------
% % Setup problem parameters
%% Exact solution 1 - 
u_exact = @(x,y) x.*(1-x).*y.*(1-y);
% gradient of exact solution
gradu_x_exact = @(x,y) (1-2*x)*y*(1-y);
gradu_y_exact = @(x,y) x*(1-x)*(1-2*y);
lap_exact=@(x,y) -2*( y-y^2 + x-x^2);
% Source term - 
d_exact = 1/45;
source = @(x,y) -(1+d_exact)*lap_exact(x,y);

tic;
maxTn=4;  k=[0:1];
H1err=zeros(maxTn,1);OC_H1=zeros(maxTn-1,1); L2err=zeros(maxTn,1); OC_L2=zeros(maxTn-1,1);
d_err=zeros(maxTn,1);OC_d=zeros(maxTn-1,1);
h_l=zeros(maxTn,1);NDOF=zeros(maxTn,1);nIters=zeros(maxTn,1);
for meshtype=1:2
for p_iter=1:length(k)
    Meshes=strings(1,maxTn);
    for level=1:maxTn
        tic; % ---------------HHO DataStructure------------------------------
        if(meshtype==1)
            mesh = strcat('mesh1_',num2str(level));  % triangular mesh; first two meshes are coarse meshes
        elseif(meshtype==2)
            mesh = strcat('mesh2_',num2str(level));     % Cartesian mesh; first two meshes are coarse meshes
        elseif(meshtype==3)
            mesh = strcat('hexa1_',num2str(level));     % hexa mesh
        else
            mesh = strcat('mesh4_1_',num2str(level));   % Kershaw mesh. Note that h are not halved.
        end
        Meshes(level)=mesh;
        K=k(p_iter);  %The polynomial degree on the cells and edges
        mesh_directory = 'matlab_meshes/'; %Change this to the directory containing the meshes
        hho = HHO(strcat(mesh_directory, mesh), K); h_l(level) = max(hho.mesh.h_size);
        % --------------------------------------------------------------
        hho = HHO(strcat(mesh_directory, mesh), K);
        % --------------------------------------------------------------
        h_size = max(hho.mesh.h_size);
        
        %% Global Assembling
        AS_local=cell(hho.mesh.ncells,1); A_local=cell(hho.mesh.ncells,1); b=cell(hho.mesh.ncells,1);
        M_L2_local=cell(hho.mesh.ncells,1); Mx1_L2_local=cell(hho.mesh.ncells,1);
        for ci = 1:hho.mesh.ncells %ci -> cell_i
            [AS_TF,A_TF,G_KT_ci] = DiffusionOperator(hho, ci);
            AS_local{ci}=AS_TF; A_local{ci}=A_TF;
            b{ci}=LoadOperator(hho, ci, source); 
            [MR_T,MRx1_T]=ReactionOperator(hho, ci);
            M_L2_local{ci} = MR_T; Mx1_L2_local{ci} = MRx1_T;
        end
        [A_global,AS_global,L2Mass_global,L2Mx1_global,b_global,Bd_dofs] = global_A_b(hho, A_local,AS_local,M_L2_local,Mx1_L2_local,b);

        %% Solve the problem
        ncell_dofs = hho.mesh.ncells * hho.elements{1}.ncell_dofs;
        nedge_dofs = hho.mesh.nedges * hho.elements{1}.nedge_dofs;
        ntotal_dofs = ncell_dofs + nedge_dofs;
        NDOF(level)=ntotal_dofs;
        idofs = setdiff(1:ntotal_dofs,Bd_dofs);
        
        %Initial Guess
        %Iuh = HHOInterpolate(hho, u_exact); d=uh'*(aG_global*uh);
        Prev_u = DiffSolver(hho, AS_local, b_global); Prev_d = Prev_u'*(A_global*Prev_u);
        %[uh,d,nIter] = SMW_Newton_solver(hho, AS_local, A_global,AS_global, b_global, Prev_u,Prev_d);
        [uh,d,nIter] = Newton_solver(idofs, A_global,AS_global, b_global, Prev_u,Prev_d);
        %% Compute errors
        u_interp = HHOInterpolate(hho, u_exact);
        u_diff = uh - u_interp;
        H1err(level)=sqrt(u_diff'*(AS_global*u_diff))/sqrt(u_interp'*(AS_global*u_interp));
        d_err(level)=abs(d_exact-d); 
        L2err(level)=L2Error(hho,M_L2_local, uh,u_exact,'relative');
        nIters(level)=nIter;
        if(level==1)
            fprintf('%s, k = %d, h = %2.4f | H1Err = %.4e |               | L2Err = %.4e |               | Iter: %d |\n', mesh, K, h_l(level), H1err(level),L2err(level),nIter);  
        else
            OC_H1(level-1)=log(H1err(level-1)/H1err(level))/log(h_l(level-1)/h_l(level));
            OC_L2(level-1)=log(L2err(level-1)/L2err(level))/log(h_l(level-1)/h_l(level));
            OC_d(level-1)=log(d_err(level-1)/d_err(level))/log(h_l(level-1)/h_l(level));     
            fprintf('%s, k = %d, h = %2.4f | H1Err = %.4e | OC_H1 = %2.3f | L2Err = %.4e | OC_L2 = %2.3f | Iter: %d |\n', mesh, K, h_l(level), H1err(level),OC_H1(level-1), L2err(level),OC_L2(level-1),nIter);
        end
        %% Plot the solution:
        plot_HHOCellAvg_Soln(hho,uh,u_exact)
    end
end
end





