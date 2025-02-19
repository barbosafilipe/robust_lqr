function [u_h,x_h,erro_h] = h_infinity(x0,N,x_ref,lim_sup,lim_inf,Ts,max_rate)

%--------------------------------------------------------------------------
% Implementation of the H_infinity controller - (Hassibi/Sayed/Kailath - pg 228-234)


% Global variables
  global n F_ref G_ref Q R H Ef Eg P0 F_inc G_inc un_rate
  
Hor = N;
PN = P0;  %-----> Ponderação estado final

gamma = 14350;

G1 = H;
G2 = G_ref;

L = eye(n);                 % s = Lx - regulação de s

Pc(:,:,1) = PN;

Q_h  = [1 0 0 0 0 0;...
        0 1 0 0 0 0;...
        0 0 1 0 0 0;...
        0 0 0 1 0 0;...
        0 0 0 0 1268 0;...
        0 0 0 0 0 10];
    
R_h  = [25660 0; 0 25660];

Qc = R;
Rc = Q;

Qw = eye(1);

for i = 1:Hor+1
    
    Rce = [Qc + G2'*Pc(:,:,i)*G2 G2'*Pc(:,:,i)*G1;...
           G1'*Pc(:,:,i)*G2 -(gamma^2)*Qw+G1'*Pc(:,:,i)*G1];
    Kc(:,:,i) = inv(Rce)*[G2';G1']*Pc(:,:,i)*F_ref;
    Pc(:,:,i+1) = F_ref'*Pc(:,:,i)*F_ref + L'*Rc*L - Kc(:,:,i)'*Rce*Kc(:,:,i);
    Pc(:,:,i+1) = 0.5*(Pc(:,:,i+1)+Pc(:,:,i+1)');
    
end
%     Pc(:,:,i+1)
%     eig(Pc(:,:,i+1))
%% Conditions of existence

x_h(:,1) = x0;
z = [0;0];
PI_o = eye(n);

if eig(inv(PI_o)-gamma^(-2)*Pc(:,:,Hor+2))>0
    
    for i = 1:Hor+1
        Delta(:,i) = -(gamma^2)*Qw + G1'*Pc(:,:,Hor+2-i)*G1 - ...
                    G1'*Pc(:,:,Hor+2-i)*G2*inv(Qc +...
                    G2'*Pc(:,:,Hor+2-i)*G2)*G2'*Pc(:,:,Hor+2-i)*G1;
    
    if Delta(:,i)<0
        
    else
        display('Condição 2 falhou')
    end
    end    
    
%% Regulation
    
%adaptation to de length
    x_ref(:,3002) = 1.0e-06*[-0.2373; 0.0053; -0.0480; 0.0117; 0.0468; 0.0009];
    x_ref(:,3003) = 1.0e-06*[-0.2373; 0.0053; -0.0480; 0.0117; 0.0468; 0.0009];


%--------------------------------------------------------------------------

    for k=1:1
        
        for i = 1:Hor+1
            erro_h(:,i) = x_h(:,i)-x_ref(:,i);
            
            R_Gc = Qc + G2'*Pc(:,:,Hor+2-i)*G2;
            w(:,i) = [Ef Eg]*[erro_h(:,i);z];
%             u_h(:,i) = -inv(R_Gc)*G2'*Pc(:,:,Hor+2-i)*F*erro_h(:,i) - ...
%                         inv(R_Gc)*G2'*Pc(:,:,Hor+2-i)*G1*w(:,i);
            u_h(:,i) = -inv(R_Gc)*G2'*Pc(:,:,Hor+2-i)*erro_h(:,i)+...
                        inv(R_Gc)*G2'*Pc(:,:,Hor+2-i)*G2*z;
            u_h(1,i) = min(max(lim_inf,u_h(1,i)),lim_sup);
            u_h(2,i) = min(max(lim_inf,u_h(2,i)),lim_sup);
            
                if(i>1 && ((u_h(1,i)-u_h(1,i-1))/Ts)>max_rate)
                    u_h(1,i) = Ts*max_rate+u_h(1,i-1);
                    u_h(2,i) = Ts*max_rate*u_h(2,i-1);
                elseif(i>1 && (((u_h(1,i)-u_h(1,i-1))/Ts)<-max_rate))
                    u_h(1,i) = Ts*(-max_rate)+u_h(1,i-1);
                    u_h(2,i) = Ts*(-max_rate)+u_h(2,i-1);
                end
            
%             x_h(:,i+1) = F_inc*x_h(:,i)+[ G_inc G1 ]*[ u_h(:,i); w(:,i)];
            x_h(:,i+1) = F_inc*x_h(:,i)+G_inc*u_h(:,i);
            z = u_h(:,i);
            w_hat(:,i) = -[0 0 1]*Kc(:,:,Hor+2-i)*x_h(:,i);

        end
        
    SN = 0;
    SD = 0;
    Sp = 0;
        
        for i = 1:Hor+1
            
           SN = SN+(u_h(:,i)'*Qc*u_h(:,i)+x_h(:,i)'*L'*Rc*L*x_h(:,i));  
           SD = SD+w(:,i)*Qw*w(:,i);
           Sp = Sp+(w(:,i)-w_hat(:,i))'*Delta(:,i)*(w(:,i)-w_hat(:,i));
           
        end
        
    Quoc = (x_h(:,Hor+2)'*Pc(:,:,1)*x_h(:,Hor+2)+SN)/(x_h(:,1)'*...
            inv(PI_o)*x_h(:,1)+SD);   
    Quoc < gamma^2;
    J(:,k) = x0'*(inv(PI_o)-gamma^(-2)*Pc(:,:,Hor+2))*x0-(gamma^(-2))*Sp;
        
    end

    
else
    
    display('Condição 1 falhou')
    
end
end