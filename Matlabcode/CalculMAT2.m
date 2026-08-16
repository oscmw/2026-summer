clear; clc;

%% Geometry
Ri = 0.01;
Ro = 0.015;

%% Material
E  = 1e6;
nu = 0.42;

mu = E/(2*(1+nu));
lambda = E*nu/((1+nu)*(1-2*nu));

%% Pressure
Pi = 1e5;
Po = 0;

%% Initial mesh
R = linspace(Ri,Ro,300);


%% Initial guess
lambda_z0 = 1.0;

solinit = bvpinit(R,@guess,lambda_z0);

%% Solve BVP
sol = bvp4c(@(R,y,p) odeCylinder(R,y,p,mu,lambda), ...
            @(ya,yb,p) bcCylinder(ya,yb,p,Pi,Po,mu,lambda,Ri,Ro), ...
            solinit);

%% Evaluate solution

Y = deval(sol,R);

r  = Y(1,:);
dr = Y(2,:);


R  = R(:)';
r  = r(:)';
dr = dr(:)';


% lambda_z is a constant parameter
lambda_z_value = sol.parameters;



%% Stretches

lambda_r = dr;

lambda_t = r./R;

lambda_z = lambda_z_value*ones(size(R));

J = lambda_r .* lambda_t .* lambda_z;

%% Cauchy stresses

sigma_r = (mu./J).*(lambda_r.^2-1) ...
        + (lambda./J).*log(J);


sigma_t = (mu./J).*(lambda_t.^2-1) ...
        + (lambda./J).*log(J);


sigma_z = (mu./J).*(lambda_z.^2-1) ...
        + (lambda./J).*log(J);



%% Displacement

u = r-R;

%% Print axial stretch

fprintf('Axial stretch lambda_z = %f\n',lambda_z(1));
%% Geting file from fem
nl = 10;               %number of layers
lt = (Ro-Ri)/nl;       %layer thickness
Positon = []
Positon(1) = Ri
for i = 2:nl+1
Positon(i) = Positon(i-1) + lt;
end
Positon
PosElem = []
PosElem(1) = Ri+lt/2
for i = 2:nl
PosElem(i) = PosElem(i-1) + lt;
end
PosElem
T = readmatrix('Psweep3afm.txt');
FEAfm = T(11,2:12);
T = readmatrix('Psweep3afd.txt');
FEAfd = T(11,2:12);
T = readmatrix('Psweep3add.txt');
FEAdd = T(11,2:12);
T = readmatrix('Psweep3adi.txt');
FEAdi = T(11,2:12);
T = readmatrix('Psweep3afi.txt');
FEAfi = T(11,2:12);
T = readmatrix('Psweep3astrm.txt');
FEAS1 = T(11,2:11);
T = readmatrix('Psweep3astryy.txt');
FEAS2 = T(11,2:11);

%% Plots

figure (1)
plot(R,u,'LineWidth',2)
hold on
plot(Positon,FEAfm,'o-','LineWidth',1)
plot(Positon,FEAfd,'x-','LineWidth',1)
plot(Positon,FEAfi,'*-','LineWidth',1)
%plot(Positon,FEAdd,'*')
%plot(Positon,FEAdi,'o')
xlabel('Reference radius R')
ylabel('Radial displacement u')
legend('MATLAB','ElemMixed','ElemDisp','ElemIncop')
grid on
hold off


figure (2)
plot(R,sigma_r,'LineWidth',2)
hold on
plot(PosElem,FEAS1,'o-')
xlabel('Reference radius R')
ylabel('Stress \sigma_r')
grid on
hold off
figure (3)
plot(R,sigma_t,'LineWidth',2)
hold on
plot(PosElem,FEAS2,'o-')
xlabel('Reference radius R')
ylabel('Stress \sigma_\theta')
grid on
hold off
figure (4)
plot(R,sigma_z,'LineWidth',2)
xlabel('Reference radius R')
ylabel('Stress \sigma_z')
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial guess
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = guess(R)

q0 = 0;

y = [R
     1
     q0];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Differential equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dydR = odeCylinder(R,y,p,mu,lambda)

r  = y(1);
rp = y(2);

% p is the constant lambda_z
lz = p;


lr = rp;
lt = r/R;


J = lr*lt*lz;


sr = mu/J*(lr^2-1) ...
     + lambda/J*log(J);


st = mu/J*(lt^2-1) ...
     + lambda/J*log(J);


sz = mu/J*(lz^2-1) ...
     + lambda/J*log(J);



%% derivative of sigma_r with respect to lambda_r

eps = 1e-8;

lr2 = lr + eps;

J2 = lr2*lt*lz;

sr2 = mu/J2*(lr2^2-1) ...
      + lambda/J2*log(J2);

dsr_dlr = (sr2-sr)/eps;


%% derivative of lambda_theta

dlt_dR = (lr-lt)/R;


%% derivative of sigma_r with respect to lambda_theta

eps2 = 1e-8;

lt2 = lt + eps2;

J3 = lr*lt2*lz;

sr3 = mu/J3*(lr^2-1) ...
      + lambda/J3*log(J3);

dsr_dlt = (sr3-sr)/eps2;


%% equilibrium equation

rhs = (lr/r)*(st-sr);

rpp = (rhs - dsr_dlt*dlt_dR)/dsr_dlr;



% axial force integral
dq = sz*r;


dydR = [rp
        rpp
        dq];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Boundary conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function res = bcCylinder(ya,yb,p,Pi,Po,mu,lambda,Ri,Ro)

lz = p;


%% Inner surface

r  = ya(1);
lr = ya(2);

lt = r/Ri;

J = lr*lt*lz;


sr_i = mu/J*(lr^2-1) ...
     + lambda/J*log(J);



%% Outer surface

r  = yb(1);
lr = yb(2);

lt = r/Ro;

J = lr*lt*lz;


sr_o = mu/J*(lr^2-1) ...
     + lambda/J*log(J);



%% Boundary conditions

res = [sr_i + Pi
       sr_o + Po
       ya(3)
       yb(3)];

end