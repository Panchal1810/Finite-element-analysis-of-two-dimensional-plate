
function [k_e, r_qe1] = Kel_plate_tri3(E,v,p,q,h, CCORD,n1)



gpLocs = [1/6,1/6; 2/3, 1/6; 1/6, 2/3];
gpWts = [1/6,1/6,1/6];

kk=zeros(6,6);
kp=zeros(6,6);
r_qe = zeros(6,1);

for i=1:length(gpWts)

s = gpLocs(i, 1); t = gpLocs(i, 2); w = gpWts(i);

%n1 for 3 noded triangle element
%n = [((1/4)*(1 - s)*(1 - t))+ ((1/4)*(1 - s)*(t + 1)), (1/4)*(s + 1)*(1 -t), (1/4)*(s + 1)*(t + 1)];
n = [(1-s-t)  s   t];


%n2 for 4 noded rectangle element 
%n = [(1/4)*(1 - s)*(1 - t), (1/4)*(s + 1)*(1 -t), (1/4)*(s + 1)*(t + 1), (1/4)*(1 - s)*(t + 1)];

%n3 for 8 noded rectangle element
%n3 = [-(1/4)*(-1+s)*(-1+t)*(1+s+t),  (1/2)*(1-(s^2))*(1-t),  -(1/4)*(1-t)*(1 -(s^2)+t+st), (1/2)*(1+s)*(1-(t^2)),  (1/4)*(1+s)*(1+t)*(-1+s+t), (1/2)*(1-(s^2))*(1+t),  (1/4)*(-1+s)*(1+s-t)*(1+t),  (1/2)*(-1+s)*(-1+(t^2))];

%dns=[((-1 + t)/4) + ((-1 - t)/4), (1 - t)/4, (1 + t)/4];
%dnt=[((-1 + s)/4) + ((1 - s)/4), (-1 - s)/4, (1 + s)/4];

dns = [-1  1  0];
dnt = [-1  0  1];


x = n*CCORD(n1,2); 
y = n*CCORD(n1,3);

%disp(CCORD(n1,2));

dxs = dns*CCORD(n1,2); dxt = dnt*CCORD(n1,2);
dys = dns*CCORD(n1,3); dyt = dnt*CCORD(n1,3);
J = [dxs, dxt; dys, dyt]; detJ = det(J);

dnx = (J(2, 2)*dns - J(2, 1)*dnt)/detJ;
dny = (-J(1, 2)*dns + J(1, 1)*dnt)/detJ;

b1 = [dnx(1,1) 0  dnx(1,2)  0   dnx(1,3) 0];
b2 = [0  dny(1,1) 0  dny(1,2)  0  dny(1,3)];
b3 = [dny(1,1)  dnx(1,1)  dny(1,2)  dnx(1,2)  dny(1,3)  dnx(1,3)];

b = [b1; b2; b3];
na = [n(1) 0 n(2) 0 n(3) 0; 0 n(1) 0 n(2) 0 n(3)];

d = (E/(1-(v^2)))*[1 v 0; v 1 0; 0 0 (1-v)/2];

kk = kk + detJ*w* b'*d*b;
kp = kp - detJ*w*p *na'*na;
r_qe = r_qe + detJ*w*na'*q;
end

k_e=h*kk;

r_qe1 = h*r_qe;
