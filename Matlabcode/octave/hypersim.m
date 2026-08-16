%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
% HyperTube_Simon
% Inflation of a Simon type hyperelastic tube
% ofarukbk

%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
% R  : Reference radial coord.
% a  : Deformed inner radius to be solved
% p_i: Inflation pressure
% lz : axial stretch ratio
% y_0: [A, B] Initial radii
% dW/dI = Cexp(I-I_0)

function out = hypersim(x, p_i, lz, y_0, inita,I_0)

A = y_0(1);   % A  : Undeformed inner radius
B = y_0(2);   % B  : Undeformed outer radius

trr_B	= @(a) trrbsim(x,p_i,lz,y_0,a,I_0);

%fsopt	= optimoptions(@fsolve,'Display','none','Algorithm','levenberg-marquardt');
a	= fsolve(trr_B,inita) %,fsopt); % deformed inner radius
b	= sqrt(a.^2+(B.^2-A.^2)/lz); % deformed outer radius
out	= [a, b]; % deformed radii
