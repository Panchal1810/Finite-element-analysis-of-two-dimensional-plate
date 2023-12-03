

clc
clear all;


CCORD16=xlsread('plate_re4_16.xlsx',1,'A1:D26'); %Node number and corresponding cordinates are loaded
NCA16=xlsread('plate_re4_16.xlsx',1,'G1:K17');   % Element number and connected nodes with the elemnts are loaded

NCA16B = xlsread('plate_re4_16.xlsx',1,'A30:E33');
NCA16R = xlsread('plate_re4_16.xlsx',1,'A34:E37');
NCA16T = xlsread('plate_re4_16.xlsx',1,'A38:E41');
NCA16L = xlsread('plate_re4_16.xlsx',1,'A42:E45');
%CCORD4S  =xlsread('plate_re4_4.xlsx',1,'A16:D27');

CCORD = CCORD16;
NCA = NCA16;


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
NELEMENTSX = NELEMENTS/4;
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

isol= [3:1:18, 23:1:38,43:1:50]; % Free degrees of freedom
constrain = [1 2,19,20,21,22,39,40,41,42]; % constrain degree of freedom
alldof = [1:1:50]; %all dof


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
         NCAS = NCA16B;
         beta= [0;0];
     case 2
         NCAS = NCA16R;
         beta= [F1n; F1t];
     case 3
         NCAS = NCA16T;
         beta= [0;0];
     case 4
         NCAS = NCA16L;
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



