
function B = BCAL_plate_re4_16(CCORD,n1)


pt=1/sqrt(3);
gpLocs = [-pt,-pt; -pt,pt; pt,-pt; pt,pt];
gpWts = [1,1,1,1];
B = zeros(3,8);

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

B = B + b;

end



