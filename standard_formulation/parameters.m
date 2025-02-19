%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%                        System data                                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ts = 0.01;            %sample time
    Ttotal = 30;          %trajectory time
    t = 0:Ts:Ttotal;      %time vector
    
    N = numel(t);         % horizonte
    parms = 1;            % parâmetro que caracteriza o valor de mi
                          % parms = 0 -> mu = 0 / parms = 1 -> mu ~= 0
    mi = 1e+9;            % parâmetro de penalidade ( > 0)
    alfa = 0.01;          % pertence ao parâmetro de minimização ( > 0)
                          % lambda = 1 + alfa

x0 = [0; 0; 0; 0; 0.3; -0.1];
  
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
%Centre of mass calculation - it comes from technical specifications of: 
%Scania P 360 LA6x2 R885 and Librelato trailer Tipper implement
% Tractor dimentions and weight
FOH = 1.455;
AD = 3.500;
BS = 1.300;
OAL = 7.030;
JA = 2.075;
WD = 2.600;
xf = FOH; %average front axle position 
xr = OAL-JA+(BS/2); %average rear axle position
mf = 5186; % tare front weight
mr = 3723; %tare rear weight
m_u = mf+mr; %total tare unloaden
x_cog_u = (xf*mf+xr*mr)/m_u; %unloaden tractor centre of mass

%Tipper dimentions and weight 20m³ - 8x4 application
A = 11.585;
B1 = 11.285;
B2 = 10.685;
C = 2.425;
D = 2.600;
E1 = 1.380;
E2 = 1.580;
F = 0.510;
G = 3.080;
H = 8.005;
I = 0.050;

m_tipper_u = 9370; %tere tipper mass

%Parameters values
m1 = m_u;
payload = 25000;
m2 = m_tipper_u+payload;
k1 = ((OAL^2+WD^2)/12)^(1/2);
k2 = ((A^2+D^2)/12)^(1/2);
J1 = m1*(k1^2);
J2 = m2*(k2^2);
a1 = x_cog_u-xf;
b1 = xr-x_cog_u;
a2 = 0.6*H;
b2 = 0.4*H;
l1 = a1+b1;
et1 = -0.29;
l1a = l1+et1;
l2 = a2+b2;
h1 = b1+et1;

v = 60/3.6;
Fz1 = m1*9.8*(b1/l1)-m2*9.8*(b2/l2)*(et1/l1);
Fz2 = m1*9.8*(a1/l1)+m2*9.8*(b2/l2)*(l1a/l1);
Fz3 = m2*9.8*(a2/l2);
c1 = 5.73*Fz1;
c2 = 5.73*Fz2;
c3 = 5.73*Fz3;
%---
c = c1+c2;
cs1 = a1*c1-b1*c2;
cq1 = (a1^2)*c1 + (b1^2)*c2;

cs2 = l1a*c1 +et1*c2;

% State space matrices
A1 = [-(c+c3)/v (c3*(h1+l2)-cs1-(m1+m2)*v^2)/v (c3*l2)/v c3 0 0;
      (c3*h1-cs1)/v (m2*h1*v^2-cq1-c3*h1*(h1+l2))/v -(c3*h1*l2)/v -c3*h1 0 0;
      (c3*l2)/v (m2*a2*v^2-c3*l2*(h1+l2))/v -(c3*l2^2)/v -c3*l2 0 0;
      0 0 1 0 0 0;
      1 0 0 0 0 v;
      0 1 0 0 0 0];

B1 = [c1; a1*c1; 0; 0; 0; 0];

M = [m1+m2 -m2*(h1+a2) -m2*a2 0 0 0;
    -m2*h1 J1+m2*h1*(h1+a2) m2*h1*a2 0 0 0;
    -m2*a2 J2+m2*a2*(h1+a2) J2+m2*a2^2 0 0 0;
    0 0 0 1 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1];

A = (M^-1)*A1;
B = (M^-1)*B1;

C = eye(6);
D = 0;

sys_c = ss(A,B,C,D);
[F_ref, G_ref, Cd, Dd] = ssdata(c2d(sys_c,Ts,'tustin'));%discret space state

G_ref = [G_ref G_ref]; %to optimize K, we give degrees of freedom

%% Modelo Incerto

un_rate = 1.34; %uncertainty rate
m2_inc = 9370+payload*(1+un_rate);
J2_inc = m2_inc*(k2^2);

Fz1_inc = m1*9.8*(b1/l1)-m2_inc*9.8*(b2/l2)*(et1/l1);
Fz2_inc = m1*9.8*(a1/l1)+m2_inc*9.8*(b2/l2)*(l1a/l1);
Fz3_inc = m2_inc*9.8*(a2/l2);
% c1_inc = 5.73*Fz1;
% c2_inc = 5.73*Fz2;
% c3_inc = 5.73*Fz3;
c1_inc = c1;
c2_inc = c2;
c3_inc = c3;
%---
c_inc = c1_inc+c2_inc;
cs1_inc = a1*c1_inc-b1*c2_inc;
cq1_inc = (a1^2)*c1_inc + (b1^2)*c2_inc;

cs2_inc = l1a*c1_inc +et1*c2_inc;

% State space matrices
A1_inc = [-(c_inc+c3_inc)/v (c3_inc*(h1+l2)-cs1_inc-(m1+m2_inc)*v^2)/v (c3_inc*l2)/v c3_inc 0 0;
      (c3_inc*h1-cs1_inc)/v (m2_inc*h1*v^2-cq1_inc-c3_inc*h1*(h1+l2))/v -(c3_inc*h1*l2)/v -c3_inc*h1 0 0;
      (c3_inc*l2)/v (m2_inc*a2*v^2-c3_inc*l2*(h1+l2))/v -(c3_inc*l2^2)/v -c3_inc*l2 0 0;
      0 0 1 0 0 0;
      1 0 0 0 0 v;
      0 1 0 0 0 0];

B1_inc = [c1_inc; a1*c1_inc; 0; 0; 0; 0];

M_inc = [m1+m2_inc -m2_inc*(h1+a2) -m2_inc*a2 0 0 0;
    -m2_inc*h1 J1+m2_inc*h1*(h1+a2) m2_inc*h1*a2 0 0 0;
    -m2_inc*a2 J2_inc+m2_inc*a2*(h1+a2) J2+m2_inc*a2^2 0 0 0;
    0 0 0 1 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1];

A_inc = (M_inc^-1)*A1_inc;
B_inc = (M_inc^-1)*B1_inc;

C_inc = eye(6);
D_inc = 0;

sys_c_inc = ss(A_inc,B_inc,C_inc,D_inc);
[F_inc, G_inc, Cd, Dd] = ssdata(c2d(sys_c_inc,Ts,'tustin'));%discret space state

G_inc = [G_inc G_inc]; %to optimize K, we give degrees of freedom

%% Reference Calculation

    t_lc = 0:Ts:10; %total time for lane change maneuvers
    
    u_ref = zeros(1,numel(t)); %creating u_ref vector
    u_ref(:,1000:1500) = 0.01*sin(0.4*pi*t_lc(1:501)); %turn left
    u_ref(:,2000:2500) = -0.01*sin(0.4*pi*t_lc(501:1001)); %turn righ

    u_ref = [u_ref/2;u_ref/2];

    x_ref(:,1)=[0;0;0;0;0;0];
%--------------------------------------------------------------------------
    psi_ref(1) = x_ref(6,1);
    posey_ref(1) = x_ref(5,1);
    posex_ref(1) = 0;

for k=1:numel(t)-1
    
    x_ref(:,k+1)=F_ref*x_ref(:,k)+G_ref*u_ref(:,k);
    
    psi_ref(k+1) = psi_ref(k)+x_ref(2,k)*Ts;
    velX_ref(k) = v - x_ref(1,k)*psi_ref(k+1);
    velY_ref(k) = x_ref(1,k) + v*psi_ref(k+1);
    posey_ref(k+1) = posey_ref(k) + velY_ref(k)*Ts;
    posex_ref(k+1) = posex_ref(k) + velX_ref(k)*Ts;
    x_ref(1:4,k) = [0;0;0;0];
end

H  = [1; 1; 1; 1; 1; 1];

Ef = [0.000068571865228 -0.000039852646457 -0.000011400741119 -0.000080033960575 0 -0.0066666666666667];

Eg = [-0.0666666666666667 -0.0066666666666667];
    
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %

% Dimensoes do problema

n = size(F_ref,1);
m = size(G_ref,2);
p = size(H,2);
l = size(Ef,1);
q = n+l;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %

% Matrizes de Ponderação
% Q  = [1 0 0 0 0 0;...
%       0 1 0 0 0 0;...
%       0 0 1 0 0 0;...
%       0 0 0 1 0 0;...
%       0 0 0 0 15054 0;...
%       0 0 0 0 0 100];
%     
% R  = [4998 0; 0 4998];

Q  = [1 0 0 0 0 0;...
      0 1 0 0 0 0;...
      0 0 1 0 0 0;...
      0 0 0 1 0 0;...
      0 0 0 0 25000 0;...
      0 0 0 0 0 100];
    
R  = [67070 0; 0 67070];

% Equação de Riccati
  P0 = eye(n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%