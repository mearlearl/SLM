%%% DIMENSIONLESS SOLUTION TECHNIQUE %%%

% Solving dt = bdxx + g
% subject to B.C. u(0,t) = ux(1,t) + zeta*u(1,t) = 0
% subject to I.C. u(x,0) = varphi(x)

clear;
clc;
pause('on');

%%%%%%%%%%%%%%%%%%%%%%
%     Parameters     %
%%%%%%%%%%%%%%%%%%%%%%

% Geometric Parameters
L = 0.000030;                      % Starting rod length in [m];
D = 0.025;                  % Diameter of part [m]
A = pi*(D)^2/4;             % Cross sectional area [m^2]

% Heat Transfer Parameters
h = 200;                     % Convective Coefficient in [W/m^2/K]

% Process Parameters
P = 140;                    % Laser Power [W]
f = 50e-6;                  % Laser Focal length [m]
V = 0.600;                  % Laser Scan Speed [m/s]
lambda0 = 1.06e-6;          % Laser Wavelength [m]
I0 = P/A;                    % Laser Intensity [W/m^2]
hatch = 86e-6;              % Hatch Spacing [m]
tl = 30e-6;                 % Layer thickness [m]
Th = 8;                      % Heating Cycle Time [s]
Tc = 24;

% Material Properties
alpha = 7.809e-6;           % Thermal Diffusivity [m^2/s]
kcond = 29.95;              % Thermal Conductivity [W/m/K]
rf = 0.582;                 % Reflectance
kappa = 3.8563*10^9;             % Extinction Coefficient
delta = 4*pi*kappa/lambda0;                 %4*pi*kappa/lambda0; % Absorption [1/m]

% Solution Parameters
Nx = 2;                   % Number of rod divisions
Nt = 50;                   % Number of time divisions
K = 20;                     % Number of waves
CycleNumber = 2000;

% Dimensionless Variables
J = Th*(1-rf)*I0*delta*alpha/kcond;
x = linspace(0,1,Nx); 
t = linspace(0,1,Nt);

%%%% Set all initial values

zeta = L*h/kcond;
betahot = alpha*Th/L^2;
betacold = alpha*Tc/L^2;

mu = muk(K,zeta);
eigenf = calleigenf(K,mu)';

measure1 = measurek(K,zeros(Nx,1),eigenf,x);
u = [];
uavg = [];

for cyc = 0:CycleNumber-1

%Superfluous on first iteration
    
zeta = L*h/kcond;
betahot = alpha*Th/L^2;
betacold = alpha*Tc/L^2;

mu = muk(K,zeta);
eigenf = calleigenf(K,mu)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This calculates the solution to the heating generated.  Varies on length 

QH = J/betahot*(4./(2*mu-sin(2*mu)))./mu./(mu.^2+(delta*L)^2).*((delta*L)*sin(mu)-mu.*cos(mu)+mu.*exp(-delta*L));
%Q2 = 5/beta*(4./(2*mu-sin(2*mu)))./mu./mu.*(-cos(mu)+cos(4/5*mu));
XH  = cell2mat(cellfun(@(ff) ff(x), eigenf,'UniformOutput',false));
XH = diag(QH)*XH;
TH = (1-exp(-betahot*t'*(mu.^2)))';
PartialsH = zeros(K,length(x),length(t));

for k = 1:K
    PartialsH(k,:,:) = XH(k,:)'*TH(k,:);
end    

uH = reshape(sum(PartialsH,1),length(x),length(t));

%calculates the cooling effect of initial condition in measure1

Q1 = (4)*(measure1.*mu)./(2*mu-sin(2*mu));
X1  = cell2mat(cellfun(@(ff) ff(x), eigenf,'UniformOutput',false));
X1 = diag(Q1)*X1;
T1 = exp(-betacold*t'*(mu.^2))';
Partials1 = zeros(K,length(x),length(t));

for k = 1:K
    Partials1(k,:,:) = X1(k,:)'*T1(k,:);
end    

u1 = reshape(sum(Partials1,1),length(x),length(t))+ uH;

%% Calculate the time average of the heating element

uavg1 = Time_Integral(u1,1/(Nt+1));

%%% Calulates the cooling starting from end of u1 (heating and cooling).

measure2 = measurek(K,u1(:,end),eigenf,x);

Q2 = (4)*(measure2.*mu)./(2*mu-sin(2*mu));
X2  = cell2mat(cellfun(@(ff) ff(x), eigenf,'UniformOutput',false));
X2 = diag(Q2)*X2;
T2 = exp(-betahot*t'*(mu.^2))';
Partials2 = zeros(K,length(x),length(t));

for k = 1:K
    Partials2(k,:,:) = X2(k,:)'*T2(k,:);
end    

u2 = reshape(sum(Partials2,1),length(x),length(t));

uavg2 = Time_Integral(u2,1/(Nt+1));

uavgnext = [(Th*uavg1+Tc*uavg2)/(Th+Tc);zeros(CycleNumber-1-Nx+2,1)];

uavg = [uavg uavgnext];

measure1 = measurek(K,u2(:,end),eigenf,x);

% Change the length

L = L + tl;

Nx=Nx+1;

x = linspace(0,1,Nx);

%u = [u u1 u2];

end

x = linspace(0,1,Nx-1);

sum = sum(uavg,2);
numb = [CycleNumber+1:-1:1]';

totalavg = sum./numb;
plot(x,totalavg)


%for j = 1:2*Nt
 %  plot(x,u(:,j),x,uavg,'r-');
   %p.LineWidth = 3;
  % axis([0 1 min(min(u)) max(max(u))+1]);
  % pause(0.05)
%end    
%pause(0.1);
%surf(u)
%plot(x,uavg,'r-')

function measures = measurek(K,phi,eigenf,x)
    % Integrates the norms and magnitudes of the Kth eigenfunction*psi
    % functions on the domain.
   
    measures = zeros(1,K);
    
    for k = 1:K
        measures(k) = integreat(phi'.*(eigenf{k}(x)),1/(length(x)+1));
    end    
   
end

function eigenf = calleigenf(K,mu)
    %%% Creates the Kth eigenfunctions and puts them in a cell

    eigenf = cell(1,K);

    for k = 1:K
           eigenf{k} =  @(x) sin(mu(k)*x);
    end    
    
end

function muks = muk(K,zeta)

    % Calculates the first K values of mu in the eigenvalues
    % Uses the newton method to find roots of a transcendental equation.
    
    f = @(x) x + zeta*tan(x);
    fp = @(x) 1 + zeta*sec(x)^2;
    muks = zeros(1,K);
    
    for i = 1:K
    
        if (zeta > 1) && (zeta < 300)
            muks(i) = newton(f,fp,0.03 + pi/2 + pi*(i-1),0.5e-8);
        end    
        if zeta >= 300
            muks(i) = pi*i;
        end 
        if zeta < 1
            muks(i) = pi/2 + pi*(i-1);
        end    
    end    
   
end