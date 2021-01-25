clc,
clear;
% close all;
% ==========================================
% =========== Initialization ===============
P1=10; % The maximum total power of user 1
P2=10; % The maximum total power of user 2
sigma2 = 1;
% ==========================================
iter=5;
RandStream.setGlobalStream(RandStream('mcg16807','seed',sum(100*clock)));
% ==========================================
% ========== Estimated Channels ============
% ==========================================
h11 = (1/sqrt(2))*(randn(1,iter)+1i*randn(1,iter));  % Secondary Channel Gain
H11 = abs(h11).^2;
 
h12 = (sqrt(1)/sqrt(2))*(randn(1,iter)+1i*randn(1,iter));  % Secondary Channel Gain
H12 = abs(h12).^2;

h21 = (1/sqrt(2))*(randn(1,iter)+1i*randn(1,iter)); 
H21 = abs(h21).^2;

h22 = (1/sqrt(2))*randn(1,iter)+1i*randn(1,iter);  
H22 = abs(h22).^2;
% ==========================================
% ==========================================

L=1000;
Pro_in=0.95;
% beta = sqrt(chi2inv(Pro_in,2)*(1/(1+1/eps1)));
alpha2=0.5;
dummy=0;
eps2=0:0.02:0.2;
for beta1=eps2(dummy+1:end)
    dummy = dummy+1
    eps1=sqrt(-beta1*log(1-Pro_in))
%     eps1 = sqrt(chi2inv(Pro_in,2)*(1/(1+1/beta1)));
    for cnt=1:iter
        shoma=0;
        h11_w = max(abs(h11(cnt))-eps1,0);%*max(1-eps1,0);
        h22_w = max(abs(h22(cnt))-eps1,0);%*max(1-eps1,0);
        H11_w = h11_w^2;%H11(cnt);%*max(1-eps1,0)^2;
        H22_w = h22_w^2;%H22(cnt);%*max(1-eps1,0)^2;
        Rth=alpha2*log2(1+P1*H11_w);
        % =========================================
        % =========================================
        
        alpha=alpha2;
        % =================================================
        % ===== Initialization - Optimization Problem =====
        % =================================================
        a1 = [abs(h11(cnt))^2; abs(h21(cnt))^2];
        a2 = [abs(h12(cnt))^2; abs(h22(cnt))^2];
        A1 = [sigma2; a1]*[sigma2; a1]';
        A2 = [sigma2; a2]*[sigma2; a2]';

        b1 = [0 abs(h21(cnt))^2]';
        b2 = [abs(h12(cnt))^2 0]';
        B1 = [sigma2; b1]*[sigma2; b1]';
        B2 = [sigma2; b2]*[sigma2; b2]';

        f1 = [h11(cnt)^2 h21(cnt)^2]'; F1 = f1*f1';
        f2 = [h12(cnt)^2 h22(cnt)^2]'; F2 = f2*f2';

        g1 = [0 h21(cnt)^2]'; G1 = g1*g1';
        g2 = [h12(cnt)^2 0]'; G2 = g2*g2';
        
        e1=[1 0]'; E1=e1*e1'; 
        e2=[0 1]'; E2=e2*e2';
        
        E_hat1 = [0 zeros(1,2); zeros(2,1) E1];
        E_hat2 = [0 zeros(1,2); zeros(2,1) E2];
        
        K1 = [0 0.5*e1'; 0.5*e1 zeros(2,2)];
        K2 = [0 0.5*e2'; 0.5*e2 zeros(2,2)];
        % ==================================================
        % ========== Improper - Non-robust =================
        % ==================================================
        if dummy==1
            Rub = log2(1+P1*H11(cnt))+log2(1+P2*H22(cnt));Rlb=0;
            while (Rub-Rlb)/Rlb>=1e-3
                R=1/2*(Rub+Rlb);

                cvx_begin sdp quiet
                   variable C(3,3) semidefinite symmetric  %toeplitz
                   variable Q(2,2) semidefinite hermitian 

                   maximize 1
                   subject to
                     trace(E_hat1*C) <= P1^2;
                     trace(E_hat2*C) <= P2^2;
                     trace(K1*C) >= 0;
                     trace(K2*C) >= 0;
                     trace(E1*Q) <= trace(E_hat1*C);
                     trace(E2*Q) <= trace(E_hat2*C);
                     C(1,1) == 1;

                     real(trace(A1*C)-trace(F1*Q)) >=sigma2^2;
                     real(trace(A2*C)-trace(F2*Q)) >=sigma2^2;

                     real(trace(B1*C)-trace(G1*Q)) >= sigma2^2;
                     real(trace(B2*C)-trace(G2*Q)) >= sigma2^2;

                     real(trace(A1*C)-trace(F1*Q)) >= ...
                       (2^(alpha*2*R))*(real(trace(B1*C)-trace(G1*Q)));

                     real(trace(A2*C)-trace(F2*Q)) >= ...
                       (2^((1-alpha)*2*R))*(real(trace(B2*C)-trace(G2*Q)));

                cvx_end

                if strcmpi(cvx_status,'solved')
                    Rlb=R;
                    CC_n = C;
                    QQ_n = Q;
                elseif strcmpi(cvx_status,'infeasible') %|| strcmpi(cvx_status,'failed')
                    Rub=R;
                else
    %                 keyboard;
                    Rlb=R;
                    CC=C;
                    QQ=Q;
                end
            end
            CC_t_n(dummy,cnt,:,:) = CC_n;
            QQ_t_n(dummy,cnt,:,:) = QQ_n;
            CC_n(:,:) = CC_t_n(dummy,cnt,:,:);
            QQ_n(:,:) = QQ_t_n(dummy,cnt,:,:);
            if rank(CC_n) == 1 %|| CC_svd_n(1)/CC_svd_n(2)>10^6
                Rt_imp(dummy,cnt)=R;
            else
                % ====== Generating Random Vectors =======
                c_l = mvnrnd([0 0 0],CC_n,L);
                beta_l = mvnrnd([0 0],QQ_n,L);
                % ========================================
                for ll=1:L                    
                    t_l = c_l(ll,1);
                    cl_n(ll,1) = max(0,min(P1, c_l(ll,2)/t_l));
                    cl_n(ll,2) = max(0,min(P2, c_l(ll,3)/t_l));
                    % =====================================
                    ql = beta_l(ll,:)/t_l;
                    q1_hat_n(ll,1) = min(1,cl_n(ll,1)/abs(ql(1)))*ql(1);
                    q1_hat_n(ll,2) = min(1,cl_n(ll,2)/abs(ql(2)))*ql(2);
                    % =====================================
                    % =====================================
                    C_hat=[1 cl_n(ll,:)]'*[1 cl_n(ll,:)];
                    Q_hat=q1_hat_n(ll,:)'*q1_hat_n(ll,:);
                    % =====================================
                    % =====================================

                    % =====================================
                    % =====================================
                    cte1  = (cl_n(ll,1)*H11(cnt) + cl_n(ll,2)*abs(h21(cnt))^2 + sigma2)^2 ...
                        - abs(q1_hat_n(ll,1)*h11(cnt)^2 + q1_hat_n(ll,2)*h21(cnt) ^2)^2;

                    cte2 = (cl_n(ll,2)*abs(h21(cnt))^2 + sigma2).^2 ...
                        - abs(q1_hat_n(ll,2)*h21(cnt) ^2).^2;

                    R1_i_n(ll) = .5 * log2(cte1/cte2);clear cte
                    % ====================================================
                    % ====================================================
                    cte1  = (cl_n(ll,2)*H22(cnt) + cl_n(ll,1)*abs(h12(cnt))^2 + sigma2)^2 ...
                        - abs(q1_hat_n(ll,2)*h22(cnt)^2 + q1_hat_n(ll,1)*h12(cnt) ^2)^2;

                    cte2 = (cl_n(ll,1)*abs(h12(cnt))^2 + sigma2).^2 ...
                        - abs(q1_hat_n(ll,1)*h12(cnt) ^2).^2;

                    R2_i_n(ll) = .5 * log2(cte1/cte2); clear cte
                    if alpha==0
                        RR_ap_n(ll)= R1_i_n(ll);
                    elseif alpha==1
                        RR_ap_n(ll)= R2_i_n(ll);                        
                    else
                        RR_ap_n(ll)= min(R1_i_n(ll)/alpha,...
                            R2_i_n(ll)/(1-alpha));                        
                    end
                end
                [Rt_imp_n(dummy), num_n] = max(RR_ap_n);
                q1_prime_n(cnt,:)=q1_hat_n(num_n,:);
                cl_p_n(cnt,:)=cl_n(num_n,:);

            end
        end
        % =====================================
        % =====================================
        h21_w = abs(h21(cnt)) + eps1;
        h12_w = abs(h12(cnt)) + eps1;
        % =================================
        theta1 = 2*asin(min(eps1/abs(h21(cnt)),1));
%         theta2 = 2*asin(min(eps1/abs(h11(cnt)),1));
%         theta1 = theta1 + theta2;
        phi_1  = -theta1:0.01:theta1;
        cte = abs(q1_prime_n(cnt,1)*h11_w^2+q1_prime_n(cnt,2)*h21_w ^2.*...
            exp(1i.*(phi_1 + 2*angle(h21(cnt))-2*angle(h11(cnt))))).^2;

        cte1  = (cl_p_n(cnt,1)*H11_w + cl_p_n(cnt,2)*abs(h21_w)^2 + sigma2).^2 ...
            - max(cte);%abs(q1_hat(ll,1)*h11^2 + q1_hat(ll,2)*h21_w ^2).^2;

        cte2 = (cl_p_n(cnt,2)*abs(h21_w)^2 + sigma2).^2 ...
            - abs(q1_prime_n(cnt,2)*h21_w ^2).^2;

        R1_imp_n(dummy,cnt) = .5 * log2(cte1/cte2);clear cte
        % ====================================================
        % ====================================================
        theta2 = 2*asin(min(eps1/abs(h12(cnt)),1));
%         theta1 = 2*asin(min(eps1/abs(h22(cnt)),1));
%         theta2 = theta1 + theta2;
        phi_2  = -theta2:0.01:theta2;

        cte = abs(q1_prime_n(cnt,2)*h22_w^2+q1_prime_n(cnt,1)*h12_w ^2.*...
            exp(1i.*(phi_2 + 2*angle(h12(cnt))-2*angle(h22(cnt))))).^2;

        cte1  = (cl_p_n(cnt,2)*H22_w + cl_p_n(cnt,1)*abs(h12_w)^2 + sigma2).^2 ...
            - max(cte);%abs(q1_hat(ll,2)*h22^2 + q1_hat(ll,1)*h12_w ^2).^2;

        cte2 = (cl_p_n(cnt,1)*abs(h12_w)^2 + sigma2).^2 ...
            - abs(q1_prime_n(cnt,1)*h12_w ^2).^2;

        R2_imp_n(dummy,cnt) = .5 * log2(cte1/cte2); clear cte
        % =================================

%             R1_imp_n(shoma,cnt) = R1_i_n(num); R2_imp_n(shoma,cnt) = R2_i_n(num);
        clear RR_ap_n R1_i_n R2_i_n
        
        % ===========================================       
        
        
        clear temp R_temp R_temp1 R_temp2

%%        
        % =========================================
        % =========================================
        for alpha=0:0.01:1
            shoma=shoma+1;
            % =============================================
            % ========= U1-Proper/U2-Max-Power ============
            % =============================================
            R=alpha*log2(1+P2*abs(h22_w)^2/sigma2);
            Gm = 2^(2*R)-1;
            a1 = (abs(h11_w))^4;
            b1 = 2*(abs(h11_w)^2)*(sigma2 + P2*(abs(h21(cnt))+eps1)^2);
            % =============================================
            % =============================================
            a2 = Gm*((abs(h21(cnt))+eps1)^4)*((abs(h12(cnt))+eps1)^4) ...
                /((abs(h22_w))^4);

            b2 = 2*(sigma2 * Gm - P2*(abs(h22_w))^2)*((abs(h21(cnt))+eps1)^4)* ...
                ((abs(h12(cnt))+eps1)^2)/((abs(h22_w))^4);

            c2 = sigma2^2 + 2 * sigma2 * P2 *((abs(h21(cnt))+eps1)^2) - ...
                2 * sigma2 *P2 *((abs(h21(cnt)) + eps1)^4)/(abs(h22_w)^2) + ...
                Gm * sigma2^2 *((abs(h21(cnt))+eps1)^4)/(abs(h22_w)^4);
            % =====================================================
            % =====================================================
            P_max  = min((P2*abs(h22_w)^2/(2^R-1)-sigma2)/(abs(h12(cnt))+eps1)^2,P1);
            P_max2 = min(2*P2*abs(h22_w)^2/Gm - sigma2,P_max);
            if (b2*a1-b1*a2)*c2>=0
                if c2>=0
                    kappa2=0;
                    p1 = P_max;
                else
                    p1 = min(P_max,P_max2);
                    if P_max>=P_max2
                        kappa2 = 1;
                    else
                        temp1  = (sigma2+p1*(abs(h12(cnt))+eps1)^2)/(P2*abs(h22_w)^2);
                        temp2  = (2^(2*R)-1)*temp1^2-2*temp1;
                        kappa2 = min(sqrt(1-temp2),1);
                        clear temp1 temp2
                    end
                end
            else
                temp1 = a1*c2 + sqrt((a1*c2)^2-b1*c2*(a1*b2-a2*b1));
                temp2 = a2*b1-a1*b2;
                p_opt = temp1/temp2;
                if c2>0
                    p1 = min(P_max,p_opt);
                    if P_max > p_opt
                        temp1  = (sigma2+p1*(abs(h12(cnt))+eps1)^2)/(P2*abs(h22_w)^2);
                        temp2  = (2^(2*R)-1)*temp1^2-2*temp1;
                        kappa2 = min(sqrt(1-temp2),1);
                        clear temp1 temp2
                        if kappa2==1
                            p1 = P_max2;
                        end
                    else
                        kappa2=0;
                        p1 = P_max;                
                    end
                else
                    if P_max2 > p_opt
                        kappa2=0;
                        p1 = P_max;
                    else
                        kappa2=1;
                        p1 = P_max2;
                    end
                end

            end
            R2_imp_1 (shoma) = 0.5*log2(1 + 2*P2*abs(h22_w)^2/(sigma2+p1*(abs(h12(cnt))+eps1)^2) ...
                +(P2*abs(h22_w)^2/(sigma2+p1*(abs(h12(cnt))+eps1)^2))^2*(1-kappa2^2));
            temp1 = p1^2*abs(h11_w)^4 + 2*p1*abs(h11_w)^2*(sigma2+P2*(abs(h21(cnt))+eps1)^2);
            temp2 = sigma2^2 + 2*sigma2*P2*(abs(h21(cnt))+eps1)^2+(1-kappa2^2)*P2^2*(abs(h21(cnt))+eps1)^4;
            R1_imp_1 (shoma) = 0.5*log2(1+temp1/temp2);
            clear temp1 temp2
            p1_1 (shoma,cnt) = p1;
            kp1_1(shoma,cnt) = 0;
            p2_1 (shoma,cnt) = P2;
            kp2_1(shoma,cnt) = kappa2;
            clear p1 kappa2
            % =============================================
            % ========= U2-Proper/U1-Max-Power ============
            % =============================================
            R = (alpha)*log2(1+P1*abs(h11_w)^2/sigma2);
            Gm = 2^(2*R)-1;
            a1 = (abs(h22_w))^4;
            b1 = 2*(abs(h22_w)^2)*(sigma2 + P1*(abs(h12(cnt))+eps1)^2);
            % =============================================
            % =============================================
            a2 = Gm*((abs(h12(cnt))+eps1)^4)*((abs(h21(cnt))+eps1)^4) ...
                /((abs(h11_w))^4);

            b2 = 2*(sigma2 * Gm - P2*(abs(h11_w))^2)*((abs(h12(cnt))+eps1)^4)* ...
                ((abs(h21(cnt))+eps1)^2)/((abs(h11_w))^4);

            c2 = sigma2^2 + 2 * sigma2 * P2 *((abs(h21(cnt))+eps1)^2) - ...
                2 * sigma2 *P2 *((abs(h21(cnt)) + eps1)^4)/(abs(h22_w)^2) + ...
                Gm * sigma2^2 *((abs(h21(cnt))+eps1)^4)/(abs(h22_w)^4);
            % =====================================================
            % =====================================================
            P_max = min((P1*abs(h11_w)^2/(2^R-1)-sigma2)/(abs(h21(cnt))+eps1)^2,P2);
            P_max2 = max(min(2*P1*abs(h11_w)^2/Gm - sigma2,P_max),0);
            if (b2*a1-b1*a2)*c2>=0
                if c2>=0
                    kappa1=0;
                    p2 = P_max;
                else
                    p2 = min(P_max,P_max2);
                    if P_max>=P_max2
                        kappa1 = 1;
                    else
                        temp1  = (sigma2+p2*(abs(h21(cnt))+eps1)^2)/(P1*abs(h11_w)^2);
                        temp2  = (2^(2*R)-1)*temp1^2-2*temp1;
                        kappa1 = min(sqrt(1-temp2),1);
                        clear temp1 temp2
                    end
                end
            else
                temp1 = a1*c2 + sqrt((a1*c2)^2-b1*c2*(a1*b2-a2*b1));
                temp2 = a2*b1-a1*b2;
                p_opt = temp1/temp2;
                if c2>0
                    p2 = min(P_max,p_opt);
                    if P_max > p_opt
                        temp1  = (sigma2+p2*(abs(h21(cnt))+eps1)^2)/(P1*abs(h11_w)^2);
                        temp2  = (2^(2*R)-1)*temp1^2-2*temp1;
                        kappa1 = min(sqrt(1-temp2),1);
                        clear temp1 temp2
                        if kappa1==1
                            p2 = P_max2;
                        end
                    else
                        kappa1=0;
                        p2 = P_max;                
                    end
                else
                    if P_max2 > p_opt
                        kappa1=0;
                        p2 = P_max;
                    else
                        kappa1=1;
                        p2 = P_max2;
                    end
                end

            end
            R1_imp_2 (shoma) = 0.5*log2(1 + 2*P1*abs(h11_w)^2/(sigma2+p2*(abs(h21(cnt))+eps1)^2) ...
                +(P1*abs(h11_w)^2/(sigma2+p2*(abs(h21(cnt))+eps1)^2))^2*(1-kappa1^2));
            temp1 = p2^2*abs(h22_w)^4 + 2*p2*abs(h22_w)^2*(sigma2+P1*(abs(h12(cnt))+eps1)^2);
            temp2 = sigma2^2 + 2*sigma2*P1*(abs(h12(cnt))+eps1)^2+(1-kappa1^2)*P1^2*(abs(h12(cnt))+eps1)^4;
            R2_imp_2 (shoma) = 0.5*log2(1+temp1/temp2);
            clear temp1 temp2   
            p1_2 (shoma,cnt) = P1;
            kp1_2(shoma,cnt) = kappa1;
            p2_2 (shoma,cnt) = p2;
            kp2_2(shoma,cnt) = 0;
            clear p2 kappa1
            % =============================================
            % =========== U1-Proper-Max-Power =============
            % =============================================
            R1=log2(1+P1*abs(h11_w)^2/sigma2);
            R=alpha*R1;
            Gm   = 2^(2*R)-1;
            Gm2  = 2^(2*R1)-1;
            beta = 1-(P1*abs(h11_w)^2/sigma2)/Gm;
            % =========================================
            % ======= Improper Thershold ==============
            % =========================================
            P_p = min((P1*H11_w/(2^R-1)-sigma2)/(abs(h21(cnt))+eps1)^2,P2);
            cte1 = P1^2*H11_w^2+2*sigma2*P1*H11_w-Gm*sigma2^2;
            cte2 = 2*(abs(h21(cnt))+eps1)^2*(sigma2*Gm-P1*H11_w);
            P_im1 = cte1/cte2;
            if alpha<=0.5
                P_im1=inf;
            end
            % ==============================================
            % ==============================================
            cte1 = abs(h12(cnt)+eps1)^2*(sigma2+P1*H11_w);
            cte2 = abs(h22_w)^2*(sigma2);

            if cte1/cte2 > beta
                p2 = min(P_im1,P1);
    %             p2 = max(p2,0);
                if P_im1 <= P1
                    kappa2 = 1;
                else
                    kap=0:0.01:1;
                    cte = sqrt(beta^2+(1-kap.^2)*(Gm2/Gm-1))-beta;
                    Q = cte*sigma2./((abs(h21(cnt))+eps1)^2*(1-kap.^2));
                    temp = find(Q>p2,1);
                    if isempty(temp)
                        kappa2 = 1;
                    else
                        kappa2 = kap(temp);
                    end
                end
            else
                kappa2 = 0;
                p2 = P_p;
            end
            % =============================================
            % =============================================
            R2_imp_3 (shoma) = 0.5*log2(1 + 2*p2*abs(h22_w)^2/(sigma2+P1*(abs(h12(cnt))+eps1)^2) ...
                +(p2*abs(h22_w)^2/(sigma2+P1*(abs(h12(cnt))+eps1)^2))^2*(1-kappa2^2));
            temp1 = P1^2*abs(h11_w)^4 + 2*P1*abs(h11_w)^2*(sigma2+p2*(abs(h21(cnt))+eps1)^2);
            temp2 = sigma2^2 + 2*sigma2*p2*(abs(h21(cnt))+eps1)^2+(1-kappa2^2)*p2^2*(abs(h21(cnt))+eps1)^4;
            R1_imp_3 (shoma) = 0.5*log2(1+temp1/temp2);
            clear temp1 temp2
            p1_3 (shoma,cnt) = P1;
            kp1_3(shoma,cnt) = 0;
            p2_3 (shoma,cnt) = p2;
            kp2_3(shoma,cnt) = kappa2;
            clear p2 kappa2
            % =============================================
            % =========== U2-Proper-Max-Power =============
            % =============================================
            R1=log2(1+P2*abs(h22_w)^2/sigma2);
            R=alpha*R1;
            Gm   = 2^(2*R)-1;
            Gm2  = 2^(2*R1)-1;
            beta = 1-(P2*abs(h22_w)^2/sigma2)/Gm;
            % =========================================
            % ======= Improper Thershold ==============
            % =========================================
            P_p  = min((P2*H22_w/(2^R-1)-sigma2)/(abs(h12(cnt))+eps1)^2,P1);
            cte1 = P2^2*H22_w^2+2*sigma2*P2*H22_w-Gm*sigma2^2;
            cte2 = 2*(abs(h12(cnt))+eps1)^2*(sigma2*Gm-P2*H22_w);
            P_im1 = cte1/cte2;
            if alpha<=0.5
                P_im1=inf;
            end
            % ==============================================
            % ==============================================
            cte1 = abs(h21(cnt)+eps1)^2*(sigma2+P2*H22_w);
            cte2 = abs(h11_w)^2*(sigma2);

            if cte1/cte2 > beta
                p1 = min(P_im1,P2);
    %             p2 = max(p2,0);
                if P_im1 <= P2
                    kappa1 = 1;
                else
                    kap=0:0.01:1;
                    cte = sqrt(beta^2+(1-kap.^2)*(Gm2/Gm-1))-beta;
                    Q = cte*sigma2./((abs(h12(cnt))+eps1)^2*(1-kap.^2));
                    temp = find(Q>p1,1);
                    if isempty(temp)
                        kappa1 = 1;
                    else
                        kappa1 = kap(temp);
                    end
                end
            else
                kappa1 = 0;
                p1 = P_p;
            end
            % =============================================
            % =============================================
            R1_imp_4 (shoma) = 0.5*log2(1 + 2*p1*abs(h11_w)^2/(sigma2+P2*(abs(h21(cnt))+eps1)^2) ...
                +(p1*abs(h11_w)^2/(sigma2+P2*(abs(h21(cnt))+eps1)^2))^2*(1-kappa1^2));
            temp1 = P2^2*abs(h22_w)^4 + 2*P2*abs(h22_w)^2*(sigma2+p1*(abs(h12(cnt))+eps1)^2);
            temp2 = sigma2^2 + 2*sigma2*p1*(abs(h12(cnt))+eps1)^2+(1-kappa1^2)*p1^2*(abs(h12(cnt))+eps1)^4;
            R2_imp_4 (shoma) = 0.5*log2(1+temp1/temp2);
            clear temp1 temp2
            p1_4 (shoma,cnt) = p1;
            kp1_4(shoma,cnt) = kappa1;
            p2_4 (shoma,cnt) = P2;
            kp2_4(shoma,cnt) = 0;
            clear p1 kappa1  
        end
        Rth2 = R1_imp_n(dummy,cnt);
        Rth_2(dummy,cnt) = Rth2;
        temp1 = find(R1_imp_1>=Rth2, 1, 'last' );
        temp2 = find(R1_imp_2>=Rth2,1);
        temp3 = find(R1_imp_3>=Rth2,1);
        temp4 = find(R1_imp_4>=Rth2,1, 'last');
        temp5 = max(R2_imp_2(temp2),R2_imp_3(temp3));
        if ~isempty(temp1)
            temp5=max(R2_imp_1(temp1),temp5);
        end
        if ~isempty(temp4)
            temp5=max(R2_imp_4(temp4),temp5);
        end
        R2_imp(dummy,cnt)=temp5;
        % =====================================
        if dummy==1
            if temp5 == R2_imp_1(temp1)
                p1_n(cnt)  = p1_1 (temp1,cnt);
                kp1_n(cnt) = kp1_1(temp1,cnt);
                p2_n(cnt)  = p2_1 (temp1,cnt);
                kp2_n(cnt) = kp2_1(temp1,cnt);
            elseif temp5 == R2_imp_2(temp2)
                p1_n(cnt)  = p1_2 (temp2,cnt);
                kp1_n(cnt) = kp1_2(temp2,cnt);
                p2_n(cnt)  = p2_2 (temp2,cnt);
                kp2_n(cnt) = kp2_2(temp2,cnt);
            elseif temp5 == R2_imp_3(temp3)
                p1_n(cnt)  = p1_3 (temp3,cnt);
                kp1_n(cnt) = kp1_3(temp3,cnt);
                p2_n(cnt)  = p2_3 (temp3,cnt);
                kp2_n(cnt) = kp2_3(temp3,cnt);
            elseif temp5 == R2_imp_4(temp4)
                p1_n(cnt)  = p1_4 (temp4,cnt);
                kp1_n(cnt) = kp1_4(temp4,cnt);
                p2_n(cnt)  = p2_4 (temp4,cnt);
                kp2_n(cnt) = kp2_4(temp4,cnt);
            else
                tt=0;
            end
            R2_imp_4n(dummy,cnt)=R2_imp(dummy,cnt);
            Rth_2_4n(dummy,cnt)=Rth_2(dummy,cnt);
        else
            % =====================================
            % =====================================
            h21_w = abs(h21(cnt)) + eps1;
            h12_w = abs(h12(cnt)) + eps1;
            % =================================
            cte1  = (p1_n(cnt)*H11_w + p2_n(cnt)*abs(h21_w)^2 + sigma2).^2 ...
                - abs(p1_n(cnt)*kp1_n(cnt)*H11_w + p2_n(cnt)*kp2_n(cnt)*h21_w ^2).^2;

            cte2 = (p2_n(cnt)*abs(h21_w)^2 + sigma2).^2 ...
                - abs(p2_n(cnt)*kp2_n(cnt)*h21_w ^2).^2;

            R2_imp_4n(dummy,cnt) = .5 * log2(cte1/cte2);clear cte
            % ====================================================
            % ====================================================
            cte1  = (p2_n(cnt)*H22_w + p1_n(cnt)*abs(h12_w)^2 + sigma2).^2 ...
                - abs(p2_n(cnt)*kp2_n(cnt)*H22_w +...
                p1_n(cnt)*kp1_n(cnt)*h12_w ^2).^2;

            cte2 = (p1_n(cnt)*abs(h12_w)^2 + sigma2).^2 ...
                - abs(p1_n(cnt)*kp1_n(cnt)*h12_w ^2).^2;

            Rth_2_4n(dummy,cnt) = .5 * log2(cte1/cte2); clear cte
            % =================================
        end
         clear temp5
        clear R1_imp_1 R1_imp_2 R1_imp_3 R1_imp_4
  %%      
        % ========================================
        % ============ Proper ====================
        % ======================================== 
        if H11_w==0 || H22_w==0
            R2_p(dummy,cnt)=0;
            R1_p(dummy,cnt)=0;
        else
            Rub=log2(1+P1*H11_w)+log2(1+P2*H22_w);Rlb=0;
            while (Rub-Rlb)/Rlb>=1e-3
                R=1/2*(Rub+Rlb);


                cvx_begin quiet
                   variables C1 C2
                   maximize 1
                   subject to
                     0 <= C1 <= P1;
                     0 <= C2 <= P2;
                     H11_w*C1 >= (sigma2+(abs(h21(cnt))+eps1)^2*C2)*(2^(alpha2*R)-1);
                     H22_w*C2 >= (sigma2+(abs(h12(cnt))+eps1)^2*C1)*(2^((1-alpha2)*R)-1);
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
        
        
    end
    R2_imp_ave(dummy) = sum(R2_imp(dummy,:))/iter;
    R2_imp_ave_4n(dummy) = sum(R2_imp_4n(dummy,:))/iter;
    R2_p_ave  (dummy) = sum(R2_p(dummy,:))/iter;
    R1_i_ave  (dummy) = sum(R1_imp_n(dummy,:))/iter;
    R2_i_ave  (dummy) = sum(R2_imp_n(dummy,:))/iter;
    Rth_1_ave (dummy) = sum(Rth_2(dummy,:))/iter;
    Rth_1_ave_4n (dummy) = sum(Rth_2_4n(dummy,:))/iter;
end
figure('DefaultAxesFontSize',40)
hold on
plot(eps2,R2_imp_ave,'LineWidth',8)
plot(eps2,R2_p_ave,'LineWidth',8)
plot(eps2,R2_imp_ave_4n,'LineWidth',8)
plot(eps2,R2_i_ave,'LineWidth',8)
grid on
xlabel('\sigma^2_e')
ylabel('Average Rate of user 2')
legend('IGS', 'PGS','N-IGS','2-IGS')
legend('Location','best')

figure('DefaultAxesFontSize',40)
hold on
plot(eps2,Rth_1_ave,'LineWidth',8)
plot(eps2,R2_p_ave,'LineWidth',8)
plot(eps2,Rth_1_ave_4n,'LineWidth',8)
plot(eps2,R1_i_ave,'LineWidth',8)
grid on
xlabel('\sigma^2_e')
ylabel('Average Rate of user 1')
legend('IGS', 'PGS','N-IGS','2-IGS')
legend('Location','best')

figure('DefaultAxesFontSize',40)
hold on
plot(eps2,Rth_1_ave+R2_imp_ave,'LineWidth',8)
plot(eps2,2*R2_p_ave,'LineWidth',8)
plot(eps2,Rth_1_ave_4n+R2_imp_ave_4n,'LineWidth',8)
plot(eps2,R1_i_ave+R2_i_ave,'LineWidth',8)
grid on
xlabel('\sigma^2_e')
ylabel('Average sum rate')
legend('IGS', 'PGS','N-IGS','2-IGS')
legend('Location','best')

figure('DefaultAxesFontSize',40)
hold on
plot(eps2,min(Rth_1_ave,R2_imp_ave),'LineWidth',8)
plot(eps2,R2_p_ave,'LineWidth',8)
plot(eps2,min(Rth_1_ave_4n,R2_imp_ave_4n),'LineWidth',8)
plot(eps2,min(R1_i_ave,R2_i_ave),'LineWidth',8)
grid on
xlabel('\sigma^2_e')
ylabel('Minimum rate of users')
legend('IGS', 'PGS','N-IGS','2-IGS')
legend('Location','best')


