
function B = BCAL_plate_re8_64(CCORD,n1)


pt=0.774597;
gpLocs = [-pt,-pt; -pt,0; -pt,pt;  0 -pt; 0  0; 0 pt; pt,-pt; pt 0; pt pt];
gpWts = [0.30862, 0.493827, 0.308642, 0.493827, 0.790123, 0.493827, 0.308642, 0.493827, 0.308642];

B = zeros(3,16);

for i=1:length(gpWts)

s = gpLocs(i, 1); t = gpLocs(i, 2); w = gpWts(i);

%n1 for 3 noded triangle element
%n1 = [((1/4)*(1 - s)*(1 - t))+ ((1/4)*(1 - s)*(t + 1)), (1/4)*(s + 1)*(1 -t), (1/4)*(s + 1)*(t + 1)];


%n2 for 4 noded rectangle element 
%n = [(1/4)*(1 - s)*(1 - t), (1/4)*(s + 1)*(1 -t), (1/4)*(s + 1)*(t + 1), (1/4)*(1 - s)*(t + 1)];

%n3 for 8 noded rectangle element
n = [-(1/4)*(-1+s)*(-1+t)*(1+s+t),  (1/2)*(-1+(s^2))*(-1+t),  (1/4)*(-1+t)*(1 -(s^2)+t+s*t), -(1/2)*(1+s)*(-1+(t^2)),  (1/4)*(1+s)*(1+t)*(-1+s+t), -(1/2)*(-1+(s^2))*(1+t),  (1/4)*(-1+s)*(1+s-t)*(1+t),  (1/2)*(-1+s)*(-1+(t^2))];



dns=[-((t-1)*(2*s+t))/4, s*(t-1), -((2*s*t)*(t-1))/4,  (1-(t^2))/2,  ((2*s+t)*(1 + t))/4,  -s*(t+1), ((2*s-t)*(1+t))/4,  ((t^2)-1)/2];

dnt=[-((s+2*t)*(s-1))/4, ((s^2)-1)/2, -((1+s)*(s-2*t))/4,  -(s+1)*t, ((1 + s)*(s+2*t))/4, (1-(s^2))/2,  ((s-1)*(s-2*t))/4,  (s-1)*t];

x = n*CCORD(n1,2); 
y = n*CCORD(n1,3);

%disp(CCORD(n1,2));

dxs = dns*CCORD(n1,2); dxt = dnt*CCORD(n1,2);
dys = dns*CCORD(n1,3); dyt = dnt*CCORD(n1,3);
J = [dxs, dxt; dys, dyt]; detJ = det(J);

dnx = (J(2, 2)*dns - J(2, 1)*dnt)/detJ;
dny = (-J(1, 2)*dns + J(1, 1)*dnt)/detJ;

b1 = [dnx(1,1) 0  dnx(1,2)  0   dnx(1,3) 0 dnx(1,4)  0  dnx(1,5)  0  dnx(1,6)  0  dnx(1,7)  0  dnx(1,8)  0];
b2 = [0  dny(1,1) 0  dny(1,2)  0  dny(1,3)   0   dny(1,4) 0  dny(1,5)  0  dny(1,6)  0  dny(1,7)  0  dny(1,8)];
b3 = [dny(1,1)  dnx(1,1)  dny(1,2)  dnx(1,2)  dny(1,3)  dnx(1,3)  dny(1,4)  dnx(1,4)  dny(1,5)  dnx(1,5)  dny(1,6)  dnx(1,6)  dny(1,7)  dnx(1,7)   dny(1,8)  dnx(1,8)];

b = [b1; b2; b3];

B = B + b;

end



