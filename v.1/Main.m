%Modelling algorithm for hexahedral elements
clear all
close all
clc
%%%%%%%%%%%%%%%%%%
%%%%%IMPORTANT%%%%
%%%%%%%%%%%%%%%%%%

% 1) To run the modelling algorithm with a certain FE model, the user should only modify the section "LOAD MODEL DATA" 
% 2) The modelling algorithm is included with Abaqus calls and Python files for data extraction, as we implemented the bone adaptation algorithm via Abaqus user subroutines. The user may replace these calls, provided that the data returned to MATLAB maintains the same format.
% 3) Remember to update directory (path) in UMAT and .py files (or equivalent calls if the user modifies them)
% 4) The Abaqus UMAT subroutine requires identification of the boundary elements. Run the first part of this matlab code (until "START ITERATIONS") to obtain the variable "Idelext" (and associated .inp file).


%%%%%%%%%%%%%%%%%%%%%%%
%%%%%MAIN VARIABLES%%%%
%%%%%%%%%%%%%%%%%%%%%%%
%nn          (integer)                 : number of total nodes in the model
%nnbound     (integer)                 : number of boundary nodes (nodes on boundary surfaces)
%nel         (integer)                 : number of elements where the growth algorithm is applied
%inde1       (integer)                 : number of boundary elements (with any number of nodes on the boundary surface)
%inde3       (integer)                 : number of elements with 4 boundary nodes
%niter       (integer)                 : number of iterations for modelling algorithm
%path        (string)                  : directory of files
%node        (matrix, nn x 4)          : all the nodes of the model (in abaqus .inp format. 1st column:node number; 2nd, 3rd and 4th columns: nodal coordinates (X,Y,Z))
%nbound      (vector, nnbound x 1)     : nodes of the boundary surfaces of the model (format: column vector of node numbers). For example in a hollow cylinder, the boundary surfaces, both external and internal surfaces, are considered here in boundary nodes. 
%elem        (matrix, nel x 9)         : elements where the growth algorithm is applied (in abaqus .inp format. 1st column:element number. Rest:node numbers). Note: The algorithm has been tested for a model where the UMAT (user material subroutine in Abaqus) is used just in a central area of elements. At the boundary conditions (displacement constraints or concentrated loads) the solution should be carefully regarded. 
%elnbound    (matrix, inde1 x 5)       : list of boundary/exterior elements (elements with at least one boundary node). 1st column: element number. Rest: node numbers of boundary nodes in this element
%order       (matrix, inde1 x 4)       : position of the boundary nodes within an individual element. Same arrangement (elements) as in elnbound
%elnbound4   (matrix, inde3 x 5)       : elements of elnbound with 4 boundary nodes (defining the surface)
%order4      (matrix, inde3 x 4)       : position of the boundary nodes within an individual element (only for the elements of elnbound4)
%volumes_raw (matrix, (inde1 x 8) x 3) : increment of volume calculated at each integration point for each boundary element in the growth algorithm (in our case, calulated with the Abaqus UMAT subroutine)
%normext     (matrix, nnbound x 3)     : outward normal at the nodes, calculated through a weighted average as explained in the paper
%disp        (matrix, nnbound x 3)     : displacement calculated at each node to achieve the desired increments of volume. This is the result of the growth algorithm
%remesh      (matrix, nn x 4)          : displacement imposed to the nodes of the model for remeshing (displacement of the mesh when "disp" are applied as bounday conditions in boundary nodes)

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%LOAD MODEL DATA%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
path_RESULTS='D:/Modelling_algorithm/Results';
path='D:/Modelling_algorithm'; %define directory
niter=1; %number of iterations

dir=[path,'\Nodes.inp'];
node=load(dir);    %load all the nodes of the model (in abaqus .inp format)
nodeorig=node;     %save original nodes, as we will remesh later
[nn,~]=size(node); %obtain number of nodes

for i=1:nn         %create a pointer to node number (just in case the node numbering is unsorted)
    pointn(node(i,1))=i;
end

%load files with boundary nodes (in external and internal surfaces)
dir=[path,'\next.inp'];
next_raw=load(dir); %external surface (periosteum, for example)

dir=[path,'\nint.inp'];
nint_raw=load(dir); %internal surface (endosteum, for example)
% 
nbound=[next_raw';nint_raw']; %boundary nodes, join outer nodes of external and internal surfaces in the same vector, as outward normal will be calculated in the same way for both
% nbound=next_raw'; %in case there are only external nodes in boundary surfaces, comment out previous line and uncomment this
[nnbound,~]=size(nbound); %obtain size of vector of outer nodes

dir=[path,'\Remod_numbers.inp'];
elem=load(dir); %load the elements of the model where the growth algorithm is applied
[nel,nperel]=size(elem); %calculate element matrix size

for i=1:nel %create a pointer to element number (just in case the element numbering is unsorted)
    pointel(elem(i,1))=i;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%ANALYZE THE MESH AND OBTAIN CONNECTIVITY. BOUNDARY ELEMENTS AND THEIR EXTERNAL NODES%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[inde3,order4,elnbound4,elnbound,order]=Connect_hexa(nnbound,nbound,elem,nel,nperel);

%create a list of boundary elements, included those which have less than 4 boundary nodes
for i=1:length(elnbound(:,1))
    Idelext(i)=elem(pointel(elnbound(i,1)),1); 
end
[Idelexsort,index]=sort(Idelext);

%save list in a file for Abaqus
dir=[path,'\Idelext.inp'];
fid=fopen(dir,'wt'); 
fprintf(fid,'%d\n',Idelext');
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%
%%%START ITERATIONS%%%
%%%%%%%%%%%%%%%%%%%%%%

for i=1:niter
    
    %%%%%%%%%%%%%%%%
    %%%RUN ABAQUS%%%
    %%%%%%%%%%%%%%%%
    
    execute=['cd ',path,' & abaqus job=Model user=BRM interactive']; %Run model (file "Model.inp") with subroutine "BRM.for" to calculate increments of volume. 
    [status1(i),out1{i}] = system(execute);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%OBTAIN RESULTS FROM ABAQUS%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    execute=['cd ',path,' & abaqus python volumes_and_state_var.py']; 
    [status2(i),out2{i}] = system(execute); %run python script "volumes_and_state_var.py" to obtain from Abaqus .odb the change in volume in surface elements and the value of the state variables at the end of the simulation. This command creates a txt file
   
    dir=[path,'\volumes.inp'];
    volumes_raw=load(dir); %read volumes
   
    str1=['cd ',path,' & ren "volumes.inp" ']; %change name to avoid overwriting, thus allowing to save volumes for each step 
    str2=sprintf('"volumes_step%i.inp"',i);
    str=[str1,str2];
    [statusren1(i),outren1{i}] = system(str);

    %Sort volume elements as in elnbound (boundary elements)
    for j=1:length(elnbound(:,1))
        volumes_raw_ord((index(j)-1)*8+1:(index(j)-1)*8+8,:)=volumes_raw((j-1)*8+1:(j-1)*8+8,:);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%CALCULATE OUTWARD NORMAL DIRECTION%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [normext,areaelem]=normals_hexa(inde3,elem,order4,elnbound4,pointel,nnbound,nbound,pointn,node); %calculate normals
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%CALCULATE DISPLACEMENTS%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [disp,Amat,volumes,Volelem]=matsys(elnbound4,nbound,node,elem,normext,volumes_raw_ord,pointn,pointel,order4); %calculate displacements of each node to achieve the desired volume increments
    

    %%%%%%%%%%%%%%%%%%%%
    %%%CALCULATE BETA%%% This section belongs to our biological law governing tissue adaptation (bone modelling/remodelling). In our implementation, it was more convenient to compute it here rather than in the ABAQUS UMAT, but it is not part of the growth algorithm itself. In other words, this section is not required if a different biological law is used together with the growth algorithm.
    %%%%%%%%%%%%%%%%%%%%

    %calculate (area of the element/volume of the element), (A_elem/V_elem), to calculate beta, to define the modelling/remodelling fraction ((see Manchado-Morales et al. Finite element-based algorithm for organ growth and adaptation with application to bone modelling)
    l=0.04; %l=40 micrometers, in mm
    beta_pre(:,1)=elnbound4(:,1); %element number
    beta_pre(:,2)=0.04*areaelem'./Volelem; %beta value

    %write beta value for each element in a file to use it in ABAQUS UMAT
    dir=[path,'\beta_pre.inp'];
    fid1=fopen(dir,'wt');
    fprintf(fid1,'%6d, %-1.14f \n',beta_pre(1:inde3-2,:)');
    fprintf(fid1,'%6d, %-1.14f',beta_pre(inde3-1,:)');
    fclose(fid1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%END OF CALCULATE BETA%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%



    %Set the nodal displacements calculated in the previous step as boundary conditions to be imposed in the remeshing step
    Boundary=zeros(3*nnbound,4);
    for j=1:nnbound
        Boundary(3*j-2:3*j,1)=node(pointn(nbound(j)),1);
        Boundary(3*j-2,2:3)=1;
        Boundary(3*j-1,2:3)=2;
        Boundary(3*j,2:3)  =3;
        Boundary(3*j-2,4)=disp(j,1);
        Boundary(3*j-1,4)=disp(j,2);
        Boundary(3*j,4)  =disp(j,3);
    end
    
    %create Abaqus *.inp file with displacement boundary conditions
    dir=[path,'\Impdisp.inp'];
    fid1=fopen(dir,'wt');
    fprintf(fid1,'%6d, %6d, %6d, %-1.14f \n',Boundary');
    fclose(fid1);
    
    %%%%%%%%%%%%
    %%%REMESH%%%
    %%%%%%%%%%%%
    
    %run the remesh simulation
    execute=['cd ',path,' & abaqus job=Tibia_remesh interactive'];
    [status3(i),out3{i}] = system(execute);
    
    %obtain displacements from all model nodes
    execute=['cd ',path, ' & abaqus python remesh.py'];
    [status4(i),out4{i}] = system(execute);
    
    dir=[path,'\remesh.inp'];
    remesh=load(dir);  %load displacements of all model nodes
       
    str1=['cd ',path,' & ren "remesh.inp" ']; %change name to avoid overwriting, thus allowing to save the remesh in each step 
    str2=sprintf('"remesh_step%i.inp"',i);
    str=[str1,str2];
    [statusren2(i),outren2{i}] = system(str);

    %rename the old file of nodal coordinates
    str1=['cd ',path,' & ren "Nodes.inp" ']; %change name to avoid overwriting, thus allowing to save nodal coordinates in each step 
    str2=sprintf('"Nodes_step%i.inp"',i);
    str=[str1,str2];
    [statusren3(i),outren3{i}] = system(str);
    
    node(:,2:4)=node(:,2:4)+remesh(:,2:4); %update nodal coordinates

    %write new file of nodes
    dir=[path,'\Nodes.inp'];
    fid2=fopen(dir,'wt');
    fprintf(fid2,'%6d, %-1.14f, %-1.14f, %-1.14f \n',node');
    fclose(fid2);
end
    %write final results in the specified directory
    if ~exist(path_RESULTS, 'dir')
    mkdir(path_RESULTS); % create folder in case it does not exist
    end
    dir=[path_RESULTS,'\Nodes.inp'];
    fid2=fopen(dir,'wt');
    fprintf(fid2,'%6d, %-1.14f, %-1.14f, %-1.14f \n',node');
    fclose(fid2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%PLOTS (not needed in the final algorithm, just for check and test)%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%plot original boundary surface and normals in nodes
% figure
% QUAD=[nodeorig(pointn(nbound(elnbound4(:,2))),1),nodeorig(pointn(nbound(elnbound4(:,3))),1),nodeorig(pointn(nbound(elnbound4(:,4))),1), nodeorig(pointn(nbound(elnbound4(:,5))),1)];
% X=nodeorig(:,2);
% Y=nodeorig(:,3);
% Z=nodeorig(:,4);
% quadmesh(QUAD,X,Y,Z)
% for i=1:nnbound
%      px=nodeorig(pointn(nbound(i)),2);
%      py=nodeorig(pointn(nbound(i)),3);
%      pz=nodeorig(pointn(nbound(i)),4);
%      point=[px,py,pz];
%      normal=normext(i,1:3);
%      hold on
%      quiver3(px,py,pz,normal(1),normal(2),normal(3),'Color','r','LineWidth',1)
% end
% 
%  %%plot comparison between oringinal and final boundary surface 
% figure
% quadmesh(QUAD,X,Y,Z)
% hold on
% QUAD=[node(pointn(nbound(elnbound4(:,2))),1),node(pointn(nbound(elnbound4(:,3))),1),node(pointn(nbound(elnbound4(:,4))),1), node(pointn(nbound(elnbound4(:,5))),1)];
% X=node(:,2);
% Y=node(:,3);
% Z=node(:,4);
% quadmesh(QUAD,X,Y,Z)

              
                 