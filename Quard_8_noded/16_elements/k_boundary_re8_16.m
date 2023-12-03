function [ka1, rb1] = k_boundary_re8_16(side, alpha, beta, h, CCORD,n1)


% [ka, rb] = BVPQuad4NBCTerm(side, alpha, beta, coord)
% Generates kalpha and rbeta when NBC is specified along a side
% side = side over which the NBC is specified
% alpha and beta = coefficients specifying the NBC
% coord = coordinates at the element ends

% Use 2 point integration. Gauss point locations and weights


pt=0.774597;
gpLocs = [-pt, 0, pt];
gpWts = [0.555556,  0.888889,  0.555556];
ka=zeros(16,16); rb=zeros(16,1);

for i=1:length(gpWts)

a = gpLocs(i); w = gpWts(i);

switch (side)

case 1
n = [(a*(a-1))/2,  -((a^2)-1), -(a*(1 + a))/2, 0, 0, 0, 0, 0];
dna = [(2*a-1)/2,  -2*a,  -(2*a+1)/2, 0, 0, 0, 0, 0];
case 2
n = [0, 0,  (2*a*(a-1))/4, -((a^2)-1),  (a*(a+1))/2, 0, 0, 0];
dna = [0, 0, (2*a-1)/2,  -2*a,  (2*a+1)/2,  0, 0,  0];
case 3
n = [0, 0, 0, 0, -(a*(-a+1))/2, -((a^2)-1),  (a*(1+a))/2,  0];
dna = [0, 0, 0,  0, (2*a-1)/2,  -2*a,  (2*a+1)/2,  0];
case 4
n = [((1 + a)/2+((a^2)-1)/2), 0, 0, 0, 0, 0, ((1 - a)/2 + ((a^2)-1)/2),  (1-(a^2))];
dna = [(1/2)+a, 0, 0, 0, 0, 0,  (-1/2)+a, -2*a];

end

nc = [n(1) 0 n(2) 0 n(3) 0 n(4) 0 n(5)  0  n(6)  0  n(7)  0  n(8)  0; 0 n(1) 0 n(2) 0 n(3) 0 n(4) 0 n(5)  0  n(6)  0  n(7)  0  n(8)];

dxa = dna*CCORD(n1,2); dya = dna*CCORD(n1,3);
Jc=sqrt(dxa^2 + dya^2);

ny = -dxa/Jc;  nx = dya/Jc;
beta1 = [nx -ny; ny nx]*beta;

ka = ka - alpha*Jc*w*nc'*nc;
rb = rb + Jc*w*nc'*beta1;

end

ka1 = h*ka;
rb1 = h*rb;
