function [inde3,order4,elnbound4,elnbound,order]=Connect_hexa(nnbound,nbound,elem,nel,nperel)

%initialize auxiliar indices
indn=ones(1,nnbound);
inde1=1;
inde3=1;

for i=1:nel %%go through all the elements
    inde2=1;
    for j=2:nperel %go through the nodes in each element
        for k=1:nnbound %go through all the boundary nodes
            if (elem(i,j)==nbound(k))
                nboundel(k,indn(k))= elem(i,1); %to which element the node belongs to
                indn(k)=indn(k)+1; %update auxiliar index
                
                if (inde2==1)
                    elnbound(inde1,inde2)=elem(i,1); %element number
                    elnbound(inde1,inde2+1)=k; %pointer to boundary node in the element
                    order(inde1,inde2)=j-1; %position of each boundary node in the list of nodes of the element
                    inde2=inde2+1;
%                     pointeext(elem(i,1))=inde2; %pointer to boundary elements from element number
                else  
                    elnbound(inde1,inde2+1)=k; %pointer to boundary node in the element
                    order(inde1,inde2)=j-1; %position of each boundary node in the list of nodes of the element
                    inde2=inde2+1;
                end
            end
        end
    end
    %warnings    
    if (inde2>5) %element with 5 or more boundary nodes
        fprintf('Warning, element %i with 5 or more boundary nodes \n',elem(i,1))
    elseif (inde2==5) %4 exterior nodes
            elnbound4(inde3,:)=elnbound(inde1,:); %save in a new matrix those elements that have a full boundary face 
            order4(inde3,:)=order(inde1,:); %save the order of the boundary nodes in the previous elements
            if (sum(order(inde1,:))~=10) && (sum(order(inde1,:))~=26) %cases in which nodes are already ordered in a circular way (superior and inferior faces of a individual element). For the 4 remaining faces of the hexahedron, to obtain this circular arrengement, we need to switch the order of the last two nodes. It is used to define trinagles in the surface to define the normals
                elnbound4(inde3,4:5)=elnbound4(inde3,5:-1:4);
                order4(inde3,3:4)=order4(inde3,4:-1:3);
            end
            inde3=inde3+1;
    elseif (inde2<5) && (inde2>1)
        fprintf('Warning, element %i with 3 or less boundary nodes \n',elem(i,1))
    end
    if (inde2>1)
        inde1=inde1+1;
    end
end


end