function plot_HHOCellAvg_Soln(hho,uh,u_exact)

    ncell = hho.mesh.ncells;
    cell_v = hho.mesh.cell_vertices;%n4e
    vertex = hho.mesh.vertices; %c4n
    n4e=zeros(ncell,size(cell_v{1},2)-1); %s4e=zeros(ncell,hho.elements{1}.nedges); 

    for ci=1:ncell
        n4e(ci,:)=cell_v{ci}(1:hho.elements{1}.nedges);
    end

    figure
    [uh_cell,uh_edge] = HHO_Cell_Edge_Ave(hho,uh); %compute the average value of u on the cell and on the edges(for plotting)
    subplot(1,2,1)
    trisurf(n4e,vertex(:,1),vertex(:,2),u_exact(vertex(:,1),vertex(:,2)));
    subplot(1,2,2)
    ShowCellValues(vertex,n4e,uh_cell)
end