function [K] = loop_tractor(v)
persistent P
global m n p q l Q R H Ef Eg Ts Ttotal mi c c1 c2 c3 h1 l2 a2 b2 m2 cs1 cq1 a1 m1 J1 J2

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

C = eye(n);
D = 0;

sys_c = ss(A,B,C,D);
[F, G, Cd, Dd] = ssdata(c2d(sys_c,Ts,'tustin'));%discret space state

G = [G G]; %to optimize K, we give degrees of freedom


if isempty(P)
    P = eye(n);
end

lambda = norm(mi*H'*H) + 0.01;

sigma = [(mi^-1)*eye(n)-(lambda^-1)*H*H' zeros(n,l);
         zeros(l,n) ((lambda^-1)*eye(l))];

F_est = [F;Ef];
G_est = [G;Eg];
I_est = [eye(n); zeros(l,n)];

V = [ zeros(n,n)  zeros(n,m)  zeros(n,n)
      zeros(m,n)  zeros(m,m)  zeros(m,n)
      zeros(n,n)  zeros(n,m)    -eye(n)
      zeros(q,n)  zeros(q,m)    F_est
        eye(n)    zeros(n,m)  zeros(n,n)
      zeros(m,n)    eye(m)    zeros(m,n)];
  
U = [ zeros(n,n) ; zeros(m,n) ; -eye(n) ; F_est ; zeros(n,n) ; zeros(m,n)];

M = [inv(P) zeros(n,m) zeros(n,n) zeros(n,q)   eye(n)   zeros(n,m)
       zeros(m,n)   inv(R)    zeros(m,n) zeros(m,q) zeros(m,n)   eye(m)
       zeros(n,n)  zeros(n,m)   inv(Q)   zeros(n,q) zeros(n,n) zeros(n,m)
       zeros(q,n)  zeros(q,m) zeros(q,n)   sigma      I_est     -G_est
         eye(n)    zeros(n,m) zeros(n,n)   I_est'   zeros(n,n) zeros(n,m)
       zeros(m,n)    eye(m)   zeros(m,n)  -G_est'   zeros(m,n) zeros(m,m)];
 
M_cal = V'*(M\U); 
   
L = zeros(n);
K = zeros(m,n);
L(:,:) = M_cal(1:n,1:n);
K(:,:) = M_cal(n+1:n+m,1:n);
P(:,:) = M_cal(n+m+1:2*n+m,1:n);
end