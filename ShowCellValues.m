function ShowCellValues(c4n,n4e,Uh_elem)
%figure
hold on
for j=1:size(n4e,1)
    trisurf(1:size(n4e,2),c4n(n4e(j,:),1),c4n(n4e(j,:),2),ones(1,size(n4e,2))*Uh_elem(j));   
end    
view(-60,50);
hold off
end
 
