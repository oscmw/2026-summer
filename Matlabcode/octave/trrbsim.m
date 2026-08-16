%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
% This function equates the radial stress at the outer wall to zero and finds the root 'a'
function trr = trrbsim(x, p_i, lz, y_0, a, I_0)

mu = x(1); % L  : Material constant of "A" in Simon's model
k  = x(2); % k  : Material constant of "k" in Simon's model

A  = y_0(1); % A  : Undeformed inner radius
B  = y_0(2); % B  : Undeformed outer radius

%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
r   = @(R) sqrt(a.^2 + (R.^2 - A^2)/lz);	% deformed radius
lr  = @(R) R./(r(R)*lz);			% radial stretch
lt  = @(R) r(R)/R;				% circumferential stretch
I   = @(R) lr(R).^2 + lt(R).^2 + lz.^2;		% 1st invariant of deformation tensor [C]
% There are three versions of Simon Mat Model. I_0 sent as input to fcn.
% I_0 = 0;		% DONT USE THIS		%1971: I_0 does not exist
% I_0 = 3;		% USE WITH PULL-BACK	%1972: ~value of I when p_i = 0
% I_0 = lz.^2 + 2/lz;	% USE WITH DIRECT	%1972: value of I when p_i = 0

%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
% old
% intgrnd = @(R) 2.*mu.* exp(k.*( lr(R).^2 + lt(R).^2 + lz.^2)) .*(lr(R).^2-lt(R).^2).*...
%	(R./(sqrt(a.^2 + (R.^2 - A.^2)./lz).*lz))./(sqrt(a.^2 + (R.^2 - A.^2)./lz));

intgrnd = @(R) exp(k.*( I(R)-I_0 )) .* (lr(R).^2 - lt(R).^2) .* lr(R)./(r(R));
q = integral(intgrnd ,A,B,'ArrayValued',true);
trr = 2*mu*q + p_i


%ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

% trr = -p(R) + 2*L*exp(k*I)*lr^2;   % Radial stress
% ttt = -p(R) + 2*L*exp(k*I)*lt^2;   % Circumferential stress

% Equation of Equilibrium in Radial Direction in Eulerian Description
% diff(trr,R)*(1/diff(f,R))+(trr-ttt)/f = 0 % ignore this

% EoE = diff(trr,R)*(1/diff(f,R))+(trr-ttt)/f == 0; % from analytical soln
% script ignore this

% delW_delI = @(R) A*exp(k*I)
% delW_delI = A*exp(k*((R/(r(R)*lz))^2 + (r(R)/R)^2 + lz^2))
% L(R) = 2*integral(delW_delI*(lr^3-lt^2*lr)/r,R,B)
% p(R) = 2*delW_delI*lr^2 + L(R) + C
% C = -p_i - ("p(A)-C") + 2*delW_delI*lr^2
% sig_rr(A) = -p(A) + 2*delW_delI*lr^2 = p_i

