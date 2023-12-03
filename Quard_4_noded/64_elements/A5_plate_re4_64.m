

clc
clear all;


CCORD64=xlsread('plate_re4_64.xlsx',1,'A1:D82'); %Node number and corresponding cordinates are loaded
NCA64=xlsread('plate_re4_64.xlsx',1,'G1:K65');   % Element number and connected nodes with the elemnts are loaded

NCA64B = xlsread('plate_re4_64.xlsx',1,'A85:E92');
NCA64R = xlsread('plate_re4_64.xlsx',1,'A93:E100');
NCA64T = xlsread('plate_re4_64.xlsx',1,'A101:E108');
NCA64L = xlsread('plate_re4_64.xlsx',1,'A109:E116');
%CCORD4S  =xlsread('plate_re4_4.xlsx',1,'A16:D27');

CCORD = CCORD64;
NCA = NCA64;


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
NELEMENTSX = NELEMENTS/8;
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

isol= [3:1:34, 39:1:70,75:1:106,111:1:142,147:1:162]; % Free degrees of freedom
constrain = [1 2,35,36,37,38,71,72,73,74,107,108,109,110,143,144,145,146]; % constrain degree of freedom
alldof = [1:1:162]; %all dof


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
         NCAS = NCA64B;
         beta= [0;0];
     case 2
         NCAS = NCA64R;
         beta= [F1n; F1t];
     case 3
         NCAS = NCA64T;
         beta= [0;0];
     case 4
         NCAS = NCA64L;
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



