% clc,
clear;
% close all;
% % ==========================================
% % =========== Initialization ===============
P1=1; % The maximum total power of user 1
P2=1; % The maximum total power of user 2
sigma2 = 1;
% % ==========================================
iter=1;%200;

h11 = 0.0498 + 0.3731i;
H11 = abs(h11).^2;

h21 = -1.7564 + 0.5943i;
H21 = abs(h21).^2;


h12 = 0.2348 + 0.3274i;
H12 = abs(h12).^2;

h22 = -0.4498 + 0.4378i;
H22 = abs(h22).^2;
% ============================
H=[H11,H12;H21,H22]


L=1000;
Pro_in=0.95;
% beta = sqrt(chi2inv(Pro_in,2)*(1/(1+1/eps1)));
alpha2=0.5;
% alpha1=0:0.01:1;
dummy=0;
eps2=0;%0:0.02:0.2;
for t=alpha2%beta1=eps2(dummy+1:end)
    dummy = dummy+1

    for cnt=1:iter
        shoma=0;eps1=0;
        h11_w = max(abs(h11(cnt))-eps1,0);%*max(1-eps1,0);
        h22_w = max(abs(h22(cnt))-eps1,0);%*max(1-eps1,0);
        H11_w = h11_w^2;%H11(cnt);%*max(1-eps1,0)^2;
        H22_w = h22_w^2;%H22(cnt);%*max(1-eps1,0)^2;
        Rth=alpha2*log2(1+P1*H11_w);
        % =========================================
        % =========================================
        


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

%             c2 = sigma2^2 + 2 * sigma2 * P2 *((abs(h21(cnt))+eps1)^2) - ...
%                 2 * sigma2 *P2 *((abs(h21(cnt)) + eps1)^4)/(abs(h22_w)^2) + ...
%                 Gm * sigma2^2 *((abs(h21(cnt))+eps1)^4)/(abs(h22_w)^4);
            c2 = sigma2^2 + 2 * sigma2 * P2 *((abs(h21(cnt))+eps1)^2) - ...
                2 * sigma2 *P2 *((abs(h21(cnt)) + eps1)^2)/(abs(h22_w)^2) + ...
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
                
                if c2>0
                    temp1 = a1*c2 + sqrt((a1*c2)^2-b1*c2*(a1*b2-a2*b1));
                    temp2 = a2*b1-a1*b2;
                    p_opt = temp1/temp2;
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
                    temp1 = a1*c2 - sqrt((a1*c2)^2-b1*c2*(a1*b2-a2*b1));
                    temp2 = a2*b1-a1*b2;
                    p_opt = temp1/temp2;
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

%             c2 = sigma2^2 + 2 * sigma2 * P2 *((abs(h21(cnt))+eps1)^2) - ...
%                 2 * sigma2 *P2 *((abs(h21(cnt)) + eps1)^4)/(abs(h22_w)^2) + ...
%                 Gm * sigma2^2 *((abs(h21(cnt))+eps1)^4)/(abs(h22_w)^4);
            c2 = sigma2^2 + 2 * sigma2 * P2 *((abs(h21(cnt))+eps1)^2) - ...
                2 * sigma2 *P2 *((abs(h21(cnt)) + eps1)^2)/(abs(h22_w)^2) + ...
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
                
                if c2>0
                    temp1 = a1*c2 + sqrt((a1*c2)^2-b1*c2*(a1*b2-a2*b1));
                    temp2 = a2*b1-a1*b2;
                    p_opt = temp1/temp2;
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
                    temp1 = a1*c2 - sqrt((a1*c2)^2-b1*c2*(a1*b2-a2*b1));
                    temp2 = a2*b1-a1*b2;
                    p_opt = temp1/temp2;
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
                p2 = min(P_im1,P2);
    %             p2 = max(p2,0);
                if P_im1 <= P2
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
                p1 = min(P_im1,P1);
    %             p2 = max(p2,0);
                if P_im1 <= P1
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
    end

end

figure%('DefaultAxesFontSize',40)
hold on
plot(R1_imp_1,R2_imp_1,'LineWidth',8)
plot(R1_imp_2,R2_imp_2,'LineWidth',8)
%figure%('DefaultAxesFontSize',40)
%hold on
plot(R1_imp_3,R2_imp_3,'LineWidth',8)
plot(R1_imp_4,R2_imp_4,'LineWidth',8)
% plot(alpha1.*Rt_p_ave,(1-alpha1).*Rt_p_ave,'LineWidth',8)
plot([0 log2(1+P1*H11)],[log2(1+P2*H22) 0])
grid on
xlabel('Rate of user 1')
ylabel('Rate of user 2')
legend('U_1-P/U2-P_{max}', 'U_2-P/U1-P_{max}',...
    'U_1-P-P_{max}','U_2-P-P_{max}','Proper')
legend('Location','best')
% % ===========================================
