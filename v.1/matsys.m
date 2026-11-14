function [disp_ext,Amat,volumes,Volelem]=matsys(elnbound4,nbound,node,elem,normext,volumes_raw_ord,pointn,pointel,order4) %calculate displacements

%PREPARED ONLY IF ALL THE EXTERIOR ELEMENTS HAVE 4 NODES (THERE ARE NOT ELEMENTS WITH 1, 2, 3, 5, 6, 7 OR 8 EXTERNAL NODES)
%IN CASE OTHER CASES WANT TO BE CONSIDERED, ELNBOUND SHOULD BE USED INSTEAD OF ELNBOUND4 AND DEFINE HOW TO MEASURE THE VOLUMES

% nodetoGP=[1,2,4,3,5,6,8,7]; %GP assigned to each node in abaqus order
Volelem=zeros(length(elnbound4(:,1)),1);
Amat=zeros(8*length(elnbound4(:,1)),length(nbound(:,1)));
for i=1:length(elnbound4(:,1)) %go throuh all the boundary elements  WITH 4 NODES!
    %take coordinates of nodes and calculate original volumes for al the GPs
    Cart(:,:)=[node(pointn(elem(pointel(elnbound4(i,1)),2)),2:4);
               node(pointn(elem(pointel(elnbound4(i,1)),3)),2:4);
               node(pointn(elem(pointel(elnbound4(i,1)),4)),2:4);
               node(pointn(elem(pointel(elnbound4(i,1)),5)),2:4);
               node(pointn(elem(pointel(elnbound4(i,1)),6)),2:4);
               node(pointn(elem(pointel(elnbound4(i,1)),7)),2:4);
               node(pointn(elem(pointel(elnbound4(i,1)),8)),2:4);
               node(pointn(elem(pointel(elnbound4(i,1)),9)),2:4);];
    for k=1:8
        VolGP(k)=GPvol(k,Cart);
        Volelem(i)=Volelem(i)+VolGP(k); %calculate total volume element
    end
    for j=2:5 %go through all boundary nodes
        pos=order4(i,j-1);
        Cart2=Cart;
        Cart2(pos,1:3)=Cart(pos,1:3)+normext(elnbound4(i,j),:);
        for k=1:8
            VolNGP(k)=GPvol(k,Cart2);
        end
        Amat(8*(i-1)+1:8*(i-1)+8,elnbound4(i,j))= VolNGP-VolGP;         
    end
    volumes(8*(i-1)+1:8*(i-1)+8)=volumes_raw_ord(8*(i-1)+1:8*(i-1)+8,3)';
    % volumes(8*(i-1)+1:8*(i-1)+8)=volumes_raw_ord(8*(i-1)+1:8*(i-1)+8,3)'.*VolGP; %in our case, the number extracted in the ABAQUS UMAT is not multiplied by the volume and it is done here (see equation 22b in Manchado-Morales et al. Finite element-based algorithm for organ growth and adaptation with application to bone modelling)
    volumenantiguo(8*(i-1)+1:8*(i-1)+8)=VolGP;
end

%magnitudes to check
incremVOl=sum(volumes);
volorig=sum(volumenantiguo);

%depending on the problem, change the expression for the pseudoinverse
delta=(transpose(Amat)*Amat)\(transpose(Amat)*volumes');  

disp_ext=delta.*normext; %calculate displacements of each node so that volume increment of each element is the one calculate in UMAT

% for the moment there is no need to set displacements as a function of all mesh nodes because they will be applied as imposed displacements for remeshing
% disp_tot=zeros(length(ntot(:,1)),3);
% 
% for i=1:length(nbound(:,1))
%     disp_tot(nbound(i,1),:)=disp_ext(i,:);
% end


end

