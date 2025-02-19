%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%                        System's model                             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x0 = [0; 0; 0; 0; 0.3; -0.1];
  
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %

%Parameters values
m1 = 8909;
payload = 25000;
m2 = 9370+payload;
k1 = 2.16;
k2 = 3.43;
J1 = m1*(k1^2);
J2 = m2*(k2^2);
v = 60/3.6; 
a1 = 1.734;
b1 = 2.415;
a2 = 4.8;
b2 = 3.2;
l1 = a1+b1;
et1 = -0.29;
l1a = l1+et1;
l2 = a2+b2;
h1 = b1+et1;

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
[F, G, Cd, Dd] = ssdata(c2d(sys_c,Ts,'tustin'));%discret space state

G = [G G]; %to optimize K, we give degrees of freedom


% Modelo Incerto

un_rate = 1.34; %uncertainty rate
m2_inc = 9370+payload*(1+un_rate);
J2_inc = m2_inc*(k2^2);

Fz1_inc = m1*9.8*(b1/l1)-m2_inc*9.8*(b2/l2)*(et1/l1);
Fz2_inc = m1*9.8*(a1/l1)+m2_inc*9.8*(b2/l2)*(l1a/l1);
Fz3_inc = m2_inc*9.8*(a2/l2);
c1_inc = 5.73*Fz1;
c2_inc = 5.73*Fz2;
c3_inc = 5.73*Fz3;
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

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %

H  = [1; 1; 1; 1; 1; 1];

Ef = [0.000068571865228 -0.000039852646457 -0.000011400741119 -0.000080033960575 0 -0.0066666666666667];

Eg = [-0.0666666666666667 -0.0066666666666667];
    
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %

% Dimensoes do problema

n = size( F,1 );
m = size( G,2 );
p = size( H,2 );
q = size( Ef,1 );

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %

% Matrizes de Ponderação
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