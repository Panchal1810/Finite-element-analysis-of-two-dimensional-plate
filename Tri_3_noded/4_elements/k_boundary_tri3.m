
function [ka1, rb1] = k_boundary_tri3(side, alpha, beta,h, CCORD,n1)

%{
beta = [1e06;0];
CCORD = [0 1 0; 0 1 1; 0 0.5 0.5];
side =1;
n1 = [1 2 3];
alpha =0;
h = 0.02;
%}


gpLocs = [0.788675,0.211325];
gpWts = [1/2,1/2];

ka=zeros(6,6); rb=zeros(6,1);

for i=1:length(gpWts)

a = gpLocs(i); w = gpWts(i);

switch (side)

case 1
n = [(1 - a), (a), 0];
dna = [-1, 1,  0];
case 2
n = [0, (1 - a), (a)];
dna = [0, -1, 1];
case 3
n = [(a), 0, (1 - a)];
dna = [1,  0, -1];
end


nc = [n(1) 0 n(2) 0 n(3) 0; 0 n(1) 0 n(2) 0 n(3)];

dxa = dna*CCORD(n1,2); dya = dna*CCORD(n1,3);

Jc=sqrt(dxa^2 + dya^2);

ny = -dxa/Jc;  nx = dya/Jc;
beta1 = [nx -ny; ny nx]*beta;

ka = ka - alpha*Jc*w*nc'*nc;
rb = rb + Jc*w*nc'*beta1;

%disp(nc'*beta1);
%disp(Jc*w*nc'*beta1);

end
%disp('///')

ka1 = h*ka;
rb1 = h*rb;
%disp(rb1);