
function [k_e, r_qe1] = Kel_plate_re4_64(E,v,p,q,h, CCORD,n1)

%{
L = 1;
b = 1;
h = 0.02;
v = 0.25;
E = 1e11;
q = 1e6;
CCORD = [0 0; 0.5 0; 1 0; 1 0.5; 0.5 0.5; 0 0.5; 0 1; 0.5 1; 1 1];
n1 = [2 3 4 5];
%}
%; 2 3 4 5; 5 4 9 8; 6 5 8 7];

%W = CCORD(n1,1); 
%Y = CCORD(n1,2); 


% [ke, rq] BVPQuad4Element(kx, ky, p, q, coord)
% Generates for a 4 node quadrilateral element for 2d BVP
% kx, ky, p, q = parameters defining the BVP
% coord coordinates at the element ends
% Use 2x2 integr~tion. Gauss point locations and weights
pt=1/sqrt(3);
gpLocs = [-pt,-pt; -pt,pt; pt,-pt; pt,pt];
gpWts = [1,1,1,1];



kk=zeros(8,8);
kp=zeros(8,8);
r_qe=zeros(8,1);

for i=1:length(gpWts)

s = gpLocs(i, 1); t = gpLocs(i, 2); w = gpWts(i);

%n1 for 3 noded triangle element
%n1 = [((1/4)*(1 - s)*(1 - t))+ ((1/4)*(1 - s)*(t + 1)), (1/4)*(s + 1)*(1 -t), (1/4)*(s + 1)*(t + 1)];


%n2 for 4 noded rectangle element 
n = [(1/4)*(1 - s)*(1 - t), (1/4)*(s + 1)*(1 -t), (1/4)*(s + 1)*(t + 1), (1/4)*(1 - s)*(t + 1)];

%n3 for 8 noded rectangle element
%n3 = [-(1/4)*(-1+s)*(-1+t)*(1+s+t),  (1/2)*(1-(s^2))*(1-t),  -(1/4)*(1-t)*(1 -(s^2)+t+st), (1/2)*(1+s)*(1-(t^2)),  (1/4)*(1+s)*(1+t)*(-1+s+t), (1/2)*(1-(s^2))*(1+t),  (1/4)*(-1+s)*(1+s-t)*(1+t),  (1/2)*(-1+s)*(-1+(t^2))];



dns=[(-1 + t)/4, (1 - t)/4, (1 + t)/4, (-1 - t)/4];

dnt=[(-1 + s)/4, (-1 - s)/4, (1 + s)/4, (1 - s)/4];

x = n*CCORD(n1,2); 
y = n*CCORD(n1,3);

%disp(CCORD(n1,2));

dxs = dns*CCORD(n1,2); dxt = dnt*CCORD(n1,2);
dys = dns*CCORD(n1,3); dyt = dnt*CCORD(n1,3);
J = [dxs, dxt; dys, dyt]; detJ = det(J);

dnx = (J(2, 2)*dns - J(2, 1)*dnt)/detJ;
dny = (-J(1, 2)*dns + J(1, 1)*dnt)/detJ;

b1 = [dnx(1,1) 0  dnx(1,2)  0   dnx(1,3) 0 dnx(1,4)  0];
b2 = [0  dny(1,1) 0  dny(1,2)  0  dny(1,3)   0   dny(1,4)];
b3 = [dny(1,1)  dnx(1,1)  dny(1,2)  dnx(1,2)  dny(1,3)  dnx(1,3)  dny(1,4)  dnx(1,4)];

b = [b1; b2; b3];
%c = [kx, 0; 0, ky];
na = [n(1) 0 n(2) 0 n(3) 0 n(4) 0; 0 n(1) 0 n(2) 0 n(3) 0 n(4)];


d = (E/(1-(v^2)))*[1 v 0; v 1 0; 0 0 (1-v)/2];

kk = kk + detJ*w* b'*d*b;
kp = kp - detJ*w*p *na'*na;
r_qe = r_qe + detJ*w*na'*q;
end

k_e=h*kk;

r_qe1 = h*r_qe;