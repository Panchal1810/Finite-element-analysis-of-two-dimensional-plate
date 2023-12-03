

function [ka1, rb1] = k_boundary_re4(side, alpha, beta, h, CCORD,n1)


% [ka, rb] = BVPQuad4NBCTerm(side, alpha, beta, coord)
% Generates kalpha and rbeta when NBC is specified along a side
% side = side over which the NBC is specified
% alpha and beta = coefficients specifying the NBC
% coord = coordinates at the element ends

% Use 2 point integration. Gauss point locations and weights


pt=-1/sqrt(3);
gpLocs = [-pt, pt];
gpWts = [1,1];
ka=zeros(8,8); rb=zeros(8,1);

for i=1:length(gpWts)

a = gpLocs(i); w = gpWts(i);

switch (side)

case 1
n = [(1 - a)/2, (1 + a)/2, 0, 0];
dna = [-1/2, 1/2, 0, 0];
case 2
n = [0, (1 - a)/2, (1 + a)/2, 0];
dna = [0, -1/2, 1/2, 0];
case 3
n = [0, 0, (1 - a)/2, (1 + a)/2];
dna = [0, 0, -1/2, 1/2];
case 4
n = [(1 + a)/2, 0, 0, (1 - a)/2];
dna = [1/2, 0, 0, -1/2];
end

nc = [n(1) 0 n(2) 0 n(3) 0 n(4) 0; 0 n(1) 0 n(2) 0 n(3) 0 n(4)];


dxa = dna*CCORD(n1,2); dya = dna*CCORD(n1,3);

Jc=sqrt(dxa^2 + dya^2);

ny = -dxa/Jc;  nx = dya/Jc;
beta1 = [nx -ny; ny nx]*beta;

%disp(beta);
%disp(beta1);

ka = ka - alpha*Jc*w*nc'*nc;
rb = rb + Jc*w*nc'*beta1;

end

ka1 = h*ka;
rb1 = h*rb;

