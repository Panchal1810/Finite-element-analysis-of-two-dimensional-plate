clc
clear all;


CCORD4=xlsread('plate_re4_4.xlsx',1,'A1:D10'); %Node number and corresponding cordinates are loaded
NCA4=xlsread('plate_re4_4.xlsx',1,'G1:K5');   % Element number and connected nodes with the elemnts are loaded

NCA4B = xlsread('plate_re4_4.xlsx',1,'A16:E17');
NCA4R = xlsread('plate_re4_4.xlsx',1,'A18:E19');
NCA4T = xlsread('plate_re4_4.xlsx',1,'A20:E21');
NCA4L = xlsread('plate_re4_4.xlsx',1,'A22:E23');
%CCORD4S  =xlsread('plate_re4_4.xlsx',1,'A16:D27');

%{
CCORD8=xlsread('plate_re4_16.xlsx',1,'A1:D6'); %Node number and corresponding cordinates are loaded
NCA8=xlsread('plate_re4_16.xlsx',1,'G1:I5');   % Element number and connected nodes with the elemnts are loaded

CCORD16=xlsread('plate_re4_64.xlsx',1,'A1:D10'); %Node number and corresponding cordinates are loaded
NCA16=xlsread('plate_re4_64.xlsx',1,'G1:I9');   % Element number and connected nodes with the elemnts are loaded


CCORD32=xlsread('plate_re4_32.xlsx',1,'A1:D14'); %Node number and corresponding cordinates are loaded
NCA32=xlsread('plate_re4_32.xlsx',1,'G1:I13');   % Element number and connected nodes with the elemnts are loaded
%}

CCORD = CCORD4;
NCA = NCA4;


%this is the position where stress to be calculated, this is user defined
%{
[x2] = [6;8]; % for 2 elements
[x4] = [0;6;8;11]; % for 4 elements
[x8] = [0;2;3;4;6;8;10;12]; % for 8 elements
[x12] = [0;1;2;3;4;5;6;7;8;9;10;12]; % for 12 elements
%}

L = 1;
b = 1;
h = 0.02;
v = 0.25;
E = 1e11;
q = [0; 0];
p = 0;
alpha = 0;
F1n = 1e06;
F1t = 0;


% The input from the user ends here


NNODES=height(CCORD);
NELEMENTS=height(NCA);  
NELEMENTSX = NELEMENTS/2;
l=L/NELEMENTS;
DOFPN=2;          % DOF Every node has 2 dofs.
nnp=DOFPN*NNODES; % Total degree of freedom



%loading condition
K=zeros(nnp,nnp); % Intializing stiffness matrix
d=zeros(nnp,1);        %Defines size of the displacemetn matrix

R_q = zeros(nnp,1); %Initializing distributed load
R_b = zeros(nnp,1);
f = zeros(nnp,1);  %reaction matrix
Stress = zeros(3,NELEMENTS);

% for 4 elements

isol= [3 4 5 6 7 8 9 10 15 16 17 18]; % Free degrees of freedom
constrain = [1 2 11 12 13 14]; % constrain degree of freedom
alldof = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18]; %all dof


%% Calculation of element stiffness matrix
for EN=1:NELEMENTS
 %EN = 1;
  %  [B]=BCAL_beam(l,EN);

   [n1] = NCA(EN,2:5);
   
   [ ke, r_qe] = Kel_plate_re4(E,v,p,q,h,CCORD,n1); % contain kk and kp, r_q

   
    % Here, we globalize the global matrix according to the respective node
    % and respective displacement
   GNN=[2*n1(1)-1  2*n1(1)  2*n1(2)-1  2*n1(2)  2*n1(3)-1  2*n1(3) 2*n1(4)-1  2*n1(4)];

   K(GNN,GNN)=K(GNN,GNN)+ke;
 
   % force matrix
   R_q(GNN) = R_q(GNN) + r_qe;
end



for side=1:4
    
        for EN=1:NELEMENTSX
     switch (side)
     
     case 1
         NCAS = NCA4B;
         beta= [0;0];
     case 2
         NCAS = NCA4R;
         beta= [F1n; F1t];
     case 3
         NCAS = NCA4T;
         beta= [0;0];
     case 4
         NCAS = NCA4L;
         beta= [0;0];
     end

        [n1] = NCAS(EN,2:5);

      
[ke_a , r_beta] = k_boundary_re4(side,alpha,beta,h,CCORD,n1); % contain ka and rbeta

GNN=[2*n1(1)-1  2*n1(1)  2*n1(2)-1  2*n1(2)  2*n1(3)-1  2*n1(3) 2*n1(4)-1  2*n1(4)];

K(GNN,GNN)=K(GNN,GNN)+ke_a;
R_b(GNN) = R_b(GNN) + r_beta;
        end

 end



%% We have got our global matrix. Now we solve the matrices to get the
% displacements at each node
R = R_q + R_b;
d(isol)=K(isol,isol)\R(isol);
fprintf('\n----------Nodal Displacements----------');
fprintf('\nNo.  x-disp     y-disp');
for i=1:NNODES
    fprintf('\n%5d  %8.3e   %8.3e',i,d(2*i-1),d(2*i));
end


fprintf('\n');



%% Reactions at nodes

f(constrain)=K(constrain,alldof)*d;


fprintf('\n');

fprintf('\n----------Reactions----------');
fprintf('\nNo.   X-Direction   Y-direction');
for i=1:NNODES
    fprintf('\n%5d  %8.3e   %8.3e',i,f(2*i-1),f(2*i));
end

fprintf('\n');



%% Now we print the strains and stresses
fprintf('\n----------Elemental stress----------');
fprintf('\n   No.            Stress');
for EN=1:NELEMENTS

      [n1] = NCA(EN,2:5);

    GNN=[2*n1(1)-1  2*n1(1)  2*n1(2)-1  2*n1(2)  2*n1(3)-1  2*n1(3) 2*n1(4)-1  2*n1(4)];

    [B]=BCAL_plate_re4(CCORD,n1);  
   
    D = (E/(1-(v^2)))*[1 v 0; v 1 0; 0 0 (1-v)/2];

   Stress(:,EN) = D*B*d(GNN,1);  

   %Q = ['sigma_x'; 'sigma_y'; 'tau_xy'];

  fprintf('\n%5d    %8.3e',EN,Stress(:,EN));
end



