% -----------------------------------------------------
% -----------------------------------------------------
% This code derives the Pareto-optimal boundary of the rate region for
% the 2-user IC with PGS.
% ----------------------------------------------------
% ----------------------------------------------------
% ----------------------------------------------------

clc,
clear;
% close all;
% % ==========================================
% % =========== Initialization ===============
P1=10; % The maximum total power of user 1
P2=10; % The maximum total power of user 2
sigma2 = 1;
% % ==========================================
iter=1;

% Channels
h11 = 0.0498 + 0.3731i;
h21 = -1.7564 + 0.5943i;
h12 = 0.2348 + 0.3274i;
h22 = -0.4498 + 0.4378i;

H11 = abs(h11).^2;
H12 = abs(h12).^2;
H21 = abs(h21).^2;
H22 = abs(h22).^2;
H=[H11,H12;H21,H22]

L=1000;
Pro_in=0.95;
alpha2=0.5;
dummy=0;
eps2=0;%0:0.02:0.2;

% ========================================
% ============ Proper ====================
% ======================================== 
cnt=1;eps1=0;dummy=0;
for alpha2=0:0.05:1
    dummy=dummy+1

    
        Rub=log2(1+P1*H11)+log2(1+P2*H22);Rlb=0;
        while (Rub-Rlb)/Rlb>=1e-3
            R=1/2*(Rub+Rlb);


            cvx_begin quiet
               variables C1 C2
               maximize 1
               subject to
                 0 <= C1 <= P1;
                 0 <= C2 <= P2;
                 H11*C1 >= (sigma2+(abs(h21(cnt))+eps1)^2*C2)*(2^(alpha2*R)-1);
                 H22*C2 >= (sigma2+(abs(h12(cnt))+eps1)^2*C1)*(2^((1-alpha2)*R)-1);
            cvx_end    
            if strcmpi(cvx_status,'solved')
                Rlb=R;
            elseif strcmpi(cvx_status,'infeasible') || strcmpi(cvx_status,'failed')
                Rub=R;
            else
%                 keyboard;
                Rlb=R;
            end
        end
        R2_p(dummy,cnt)=(1-alpha2)*R;
        R1_p(dummy,cnt)=(alpha2)*R;clear alpha
end

plot(R1_p,R2_p,'LineWidth',2)