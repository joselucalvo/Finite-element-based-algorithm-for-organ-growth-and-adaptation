function [normext,areaelem]=normals_hexa(inde3,elem,order4,elnbound4,pointel,nnbound,nbound,pointn,node)


%initialize auxiliar indices
indn=ones(1,nnbound);

%go through all the boundary elements with 4 nodes to calculate the normal
for i=1:inde3-1
    face=order4(i,1)*1000+order4(i,2)*100+order4(i,3)*10+order4(i,4);
    switch face %check which is the boundary face (with the circular order of the boundary nodes) to determine the node in the opposite face to check the direction of the normal
        case 1234
            opface(i,:)=[5,6,7,8];
        case 1485
            opface(i,:)=[2,3,7,6];
        case 2376
            opface(i,:)=[1,4,8,5];
        case 3487
            opface(i,:)=[2,1,5,6];
        case 1265
            opface(i,:)=[4,3,7,8];
        case 5678
            opface(i,:)=[1,2,3,4];
    end
    for j=2:5 %go through the four boundary nodes
        P1=node(pointn(nbound(elnbound4(i,j))),2:4); %coordinates of the node where the normals is going to be calculated
        switch j
            case 2
                P2=node(pointn(nbound(elnbound4(i,5))),2:4); %contiguous node (node before P1 in the circular order) in this face
                P3=node(pointn(nbound(elnbound4(i,3))),2:4); %contiguous node (node after P1 in the circular order) in this face
                P4=node(pointn(nbound(elnbound4(i,4))),2:4); %node in front of P1 in this face
            case 3
                P2=node(pointn(nbound(elnbound4(i,2))),2:4);
                P3=node(pointn(nbound(elnbound4(i,4))),2:4);
                P4=node(pointn(nbound(elnbound4(i,5))),2:4);  
            case 4
                P2=node(pointn(nbound(elnbound4(i,3))),2:4);
                P3=node(pointn(nbound(elnbound4(i,5))),2:4);
                P4=node(pointn(nbound(elnbound4(i,2))),2:4); 
            case 5
                P2=node(pointn(nbound(elnbound4(i,4))),2:4);
                P3=node(pointn(nbound(elnbound4(i,2))),2:4);
                P4=node(pointn(nbound(elnbound4(i,3))),2:4); 
        end
        P5=node(pointn(elem(pointel(elnbound4(i,1)),opface(i,j-1)+1)),2:4); %node of the opposite face, in the same position as P1 in the circular order
        %calculate vectors to define area and normal
        vec1=P3-P1;
        vec2=P2-P1;
        vec3=P3-P4;
        vec4=P2-P4;
        vec5=P5-P1;
        prod=cross(vec1,vec2);
        prod2=cross(vec3,vec4);
        area(elnbound4(i,j),indn(elnbound4(i,j)))=0.5*(norm(prod)+norm(prod2)); %save the area corresponding to each node and element
        if j==2
            areaelem(i)=area(elnbound4(i,j),indn(elnbound4(i,j))); %save element area
        end
        nor=prod/norm(prod);% calculate normal
        aux=nor*vec5'; %check direction of normal
        if aux>0
            nor=-nor;
        end
        normglob(elnbound4(i,j),indn(elnbound4(i,j)),:)=nor; %save normal to access it afterwards
        indn(elnbound4(i,j))=indn(elnbound4(i,j))+1;
    end
end
[~,jx,kx]=size(normglob);
for i=1:nnbound
    Aux=zeros(jx,kx);
    Aux(:,:)=normglob(i,:,:);
    vec=zeros(1,3);
    %calculate the normal at the node as a weighted average with area 
    for j=1:(indn(i)-1)
        vec=vec+area(i,j)*Aux(j,:);
    end
    normext(i,:)=vec/norm(vec);
end


end