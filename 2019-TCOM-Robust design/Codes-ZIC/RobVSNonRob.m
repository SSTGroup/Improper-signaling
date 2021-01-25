clc,
clear;
% close all;
RandStream.setGlobalStream(RandStream('mcg16807','seed',sum(100*clock)));

% ==========================================
% =========== Initialization ===============
% P1=10; % The maximum total power of user 1
% P2=10; % The maximum total power of user 2
iter=1000; %Number of iterations

% ==========================================
% ========== Estimated Channels ============
% ==========================================
f = (1/sqrt(2))*(randn(1,iter)+1i*randn(1,iter));  
F = abs(f).^2;

g = (1/sqrt(2))*(randn(1,iter)+1i*randn(1,iter));  
G = abs(g).^2;


h = (1/sqrt(2))*(randn(1,iter)+1i*randn(1,iter));  
H = abs(h).^2;

A12 = G./F;
a12 = abs(A12);
% =====================================================
% =====================================================


bsho=0;
Asho2=0;


beta=0;
Pro_in=0.95;
Alpha=0.6;


eps2=0:.1:1;
eps2=10.^eps2;
eps1=0.02;
for P1=eps2%eps1=eps2
    P2=P1;
        
    Asho2   = Asho2+1; %counter
    
    R1_impth(Asho2)=0;
    R1_pth(Asho2)=0;
    R1_impth_non(Asho2)=0;
    R1_pth_non(Asho2)=0;
    
    iter_p(Asho2)=0;
    iter_imp(Asho2)=0;
    iter_p2(Asho2)=0;
    iter_imp2(Asho2)=0;
    out_imp_non1(Asho2)=0;
    out_p_non1(Asho2)=0;
    out_imp_non2(Asho2)=0;
    out_p_non2(Asho2)=0;
    

    beta = sqrt(-eps1*log(1-Pro_in));%sqrt(chi2inv(Pro_in,2)*(1/(1+1/eps1))); 
    
    for cnt=1:iter  
        bsho = bsho+1;
        
        
        
        % =============================================
        % ======= Robust IGS - Suboptimal =============
        % =============================================
        d_ch   = pi/2-10^(-6.5);
%         a12(cnt)    = a12(cnt)*(1+beta)^2/(1-beta)^2;
        P1ce = max(P1*H(cnt)*(1-beta)^2,0);
%         P2ce = max(P2/(G(cnt)*(1+beta)^2),0);
        % ==============================================
        % ==============================================
        R = Alpha * log2( 1 + P1ce);
        % ==============================================
        % ==============================================
        gama_R  = 2^R     - 1;
        gama_2R = 2^(2*R) - 1;
        % ==============================================
        % ==============================================
        mou=(gama_2R-P1ce-(P1ce/gama_R-1)*cos(d_ch)^2)/(gama_2R+1-sin(d_ch)^2);
        if a12(cnt)*(1+beta)^2/(1-beta)^2 > mou
            Kapa2 = 1;
            dummy = (P1ce-gama_2R)^2+((gama_2R+1)*(1-Kapa2.^2)-1+(sin(d_ch)^2)...
                *Kapa2.^2)*((1+P1ce)^2-gama_2R-1);
            dummy2 = (G(cnt)*(1+beta)^2)*...
                (-1+(gama_2R+1)*(1-Kapa2.^2)+(sin(d_ch)*Kapa2).^2);

            q2 = (P1ce-gama_2R + sqrt(dummy))./(dummy2);
            clear Kapa2
            if imag(q2)~=0 || isnan(q2) || real(q2)<0
                kapa1 = 1;
                p2 = P2;
                Kapa  = 0:0.001:0.999;
                dum  = (G(cnt)*(1+beta)^2)*((gama_2R)*(1-Kapa.^2));
                dum2 = real(((P1ce*(1+cos(d_ch).*Kapa)).^2-...
                           2*P1ce.*Kapa.*(cos(d_ch)+Kapa).*gama_2R +(gama_2R*Kapa).^2));

                Q1 = (P1ce*(1+cos(d_ch).*Kapa) - gama_2R + sqrt(dum2))./...
                           (dum);

                kapa2 = Kapa(find(Q1>=p2, 1 ));               
                if isempty(kapa2)
                    kapa2=1;
                end
                kapa1 = min( p2 * kapa2 * (G(cnt)*(1+beta)^2)/P1ce,1);
                if kapa1<1
                    Kapa2  = 0:0.001:1;

                    dummy = (P1ce-gama_2R)^2+((gama_2R+1)*(1-Kapa2.^2)-1+(sin(d_ch)^2)...
                        *Kapa2.^2)*((1+P1ce)^2-gama_2R-1);
                    dummy2 = (G(cnt)*(1+beta)^2)*...
                        (-1+(gama_2R+1)*(1-Kapa2.^2)+(sin(d_ch)*Kapa2).^2);

                    Q1 = (P1ce-gama_2R + sqrt(dummy))./(dummy2);
                    kapa2 = Kapa2(find(Q1>=P2, 1 ));
                    p2=P2;
                    if isempty(kapa2)
                        kapa2=1;
                        p2=Q1(end);
                    end                    
                end                    
            elseif q2 <= P2
                kapa2 = 1;
                p2=q2;
            else
                Kapa  = 0:0.001:1;
                Q1 = (P1ce-gama_2R+sqrt((gama_2R+1)*(P1ce^2*(1-Kapa.^2)...
                    +(gama_2R-2*P1ce)*Kapa.^2)))./((G(cnt)*(1+beta)^2)*(-1+(gama_2R+1)*(1-Kapa.^2)));
                kapa2 = Kapa(find(Q1>=P2, 1 ));
                p2=P2;
                if isempty(kapa2)
                        kapa2=1;
                        p2=Q1(end);
                end  
                
            end

        else
            kapa2=0;
            p2=min((P1ce-gama_2R+sqrt((gama_2R+1)*(P1ce^2*(1-kapa2^2)...
                +(gama_2R-2*P1ce)*kapa2^2)))/((G(cnt)*(1+beta)^2)*(-1+(gama_2R+1)*(1-kapa2^2))),P2);
            
        end

        p2  =  max(0,p2);
        kapa1 = min( cos(d_ch) * p2 * kapa2 * (G(cnt)*(1+beta)^2)/P1ce,1);
        clear Kapa2
        Kapa2 = kapa2;
        p2R=p2;kapa2R=kapa2;kappa1R=kapa1;
        %===============================================================
        %===============================================================
        %===============================================================

        
        R2_impth_2IC(Asho2,bsho)=.5 * log2(1+(p2*F(cnt)*(1-beta)^2)*((p2*F(cnt)*(1-beta)^2)*(1-Kapa2^2)+2));
        
        % =============================================
        % ============= Robust IGS ====================
        % =============================================
        theta1 = asin(beta/abs(h(cnt)));
        theta2 = asin(beta/abs(g(cnt)));
        d_ch1  = min(2*(theta1+theta2),pi);
        d_ch   = min(d_ch1,pi/2-10^(-6.5));
%         a12(cnt)    = a12(cnt)*(1+beta)^2/(1-beta)^2;
        P1ce = max(P1*H(cnt)*(1-beta)^2,0);
%         P2ce = max(P2/(G(cnt)*(1+beta)^2),0);
        % ==============================================
        % ==============================================
        R = Alpha * log2( 1 + P1ce);
        % ==============================================
        % ==============================================
        gama_R  = 2^R     - 1;
        gama_2R = 2^(2*R) - 1;
        % ==============================================
        % ==============================================
        mou=(gama_2R-P1ce-(P1ce/gama_R-1)*cos(d_ch)^2)/(gama_2R+1-sin(d_ch)^2);
        if a12(cnt)*(1+beta)^2/(1-beta)^2 > mou
            Kapa2 = 1;
            dummy = (P1ce-gama_2R)^2+((gama_2R+1)*(1-Kapa2.^2)-1+(sin(d_ch)^2)...
                *Kapa2.^2)*((1+P1ce)^2-gama_2R-1);
            dummy2 = (G(cnt)*(1+beta)^2)*...
                (-1+(gama_2R+1)*(1-Kapa2.^2)+(sin(d_ch)*Kapa2).^2);

            q2 = (P1ce-gama_2R + sqrt(dummy))./(dummy2);
            clear Kapa2
            if imag(q2)~=0 || isnan(q2) || real(q2)<0
                kapa1 = 1;
                p2 = P2;
                Kapa  = 0:0.001:0.999;
                dum  = (G(cnt)*(1+beta)^2)*((gama_2R)*(1-Kapa.^2));
                dum2 = real(((P1ce*(1+cos(d_ch).*Kapa)).^2-...
                           2*P1ce.*Kapa.*(cos(d_ch)+Kapa).*gama_2R +(gama_2R*Kapa).^2));

                Q1 = (P1ce*(1+cos(d_ch).*Kapa) - gama_2R + sqrt(dum2))./...
                           (dum);

                kapa2 = Kapa(find(Q1>=p2, 1 ));               
                if isempty(kapa2)
                    kapa2=1;
                end
                kapa1 = min( p2 * kapa2 * (G(cnt)*(1+beta)^2)/P1ce,1);
                if kapa1<1
                    Kapa2  = 0:0.001:1;

                    dummy = (P1ce-gama_2R)^2+((gama_2R+1)*(1-Kapa2.^2)-1+(sin(d_ch)^2)...
                        *Kapa2.^2)*((1+P1ce)^2-gama_2R-1);
                    dummy2 = (G(cnt)*(1+beta)^2)*...
                        (-1+(gama_2R+1)*(1-Kapa2.^2)+(sin(d_ch)*Kapa2).^2);

                    Q1 = (P1ce-gama_2R + sqrt(dummy))./(dummy2);
                    kapa2 = Kapa2(find(Q1>=P2, 1 ));
                    p2=P2;
                    if isempty(kapa2)
                        kapa2=1;
                        p2=Q1(end);
                    end                    
                end                    
            elseif q2 <= P2
                kapa2 = 1;
                p2=q2;
            else
                Kapa  = 0:0.001:1;
                Q1 = (P1ce-gama_2R+sqrt((gama_2R+1)*(P1ce^2*(1-Kapa.^2)...
                    +(gama_2R-2*P1ce)*Kapa.^2)))./((G(cnt)*(1+beta)^2)*(-1+(gama_2R+1)*(1-Kapa.^2)));
                kapa2 = Kapa(find(Q1>=P2, 1 ));
                p2=P2;
                if isempty(kapa2)
                        kapa2=1;
                        p2=Q1(end);
                end  
                
            end

        else
            kapa2=0;
            p2=min((P1ce-gama_2R+sqrt((gama_2R+1)*(P1ce^2*(1-kapa2^2)...
                +(gama_2R-2*P1ce)*kapa2^2)))/((G(cnt)*(1+beta)^2)*(-1+(gama_2R+1)*(1-kapa2^2))),P2);
            
        end

        p2  =  max(0,p2);
        kapa1 = min( cos(d_ch) * p2 * kapa2 * (G(cnt)*(1+beta)^2)/P1ce,1);
        clear Kapa2
        Kapa2 = kapa2;
        p2R=p2;kapa2R=kapa2;kappa1R=kapa1;
        %===============================================================
        %===============================================================
        %===============================================================

        
        R2_impth(Asho2,bsho)=.5 * log2(1+(p2*F(cnt)*(1-beta)^2)*((p2*F(cnt)*(1-beta)^2)*(1-Kapa2^2)+2));
        
        %===================================
        %============ Proper ===============
        %===================================
        phi2_a = 0; phi1_a = pi;
        % ===========================================
        % ======== Worst-case channels ==============
        % ===========================================
        phi2 = phi2_a - d_ch;
        phi1 = phi1_a;

        Powe_Pro=min(P2,(P1ce/gama_R-1)/(G(cnt)*(1+beta)^2));
        R2_pth(Asho2,bsho) = log2(1+Powe_Pro*F(cnt)*(1-beta)^2);
        
        cte  = (p2 * G(cnt)*(1+beta)^2 + 1)^2 - abs(p2...
            *kapa2*G(cnt)*(1+beta)^2)^2;
        
        R1_imp_th(Asho2,bsho) = .5 * log2(((p2*G(cnt)*(1+beta)^2 + P1ce + 1)^2)/(cte)...
            -(abs(p2*exp(1i*phi2)*kapa2*G(cnt)*(1+beta)^2 + P1ce*exp(1i*(phi2+pi))*kapa1).^2)/(cte));
        
        R1_p_th(Asho2,bsho)=log2(1+P1ce/(1+Powe_Pro*G(cnt)*(1+beta)^2));
        
        %%
        
        %=================================================
        %========== Non-robust Design ====================
        %=================================================
        P1ce = P1*H(cnt);
        P2ce = P2*F(cnt);
        % ==============================================
        % ==============================================
        Rnon = Alpha * log2( 1 + P1ce);
        % ==============================================
        % ==============================================
        gama_R  = 2^Rnon     - 1;
        gama_2R = 2^(2*Rnon) - 1;
        % ==============================================
        % ==============================================
        
        %mou = 1 - P1/(gama_2R - gama_R);
        if 2 * P1ce < gama_2R
            l = ((P1ce - gama_R)^2)/((gama_2R - gama_R)*(sqrt(gama_2R - 2*P1ce)-gama_R)^2);
        else
            l = 0;
        end
        % ==============================================
        % ==============================================
        aP = (P1ce/gama_R - 1)/P2ce;
        aI = (2*P1ce/gama_2R - 1)/P2ce;
        eta=.5*( aP - P2ce*aI+sqrt((aP - P2ce * aI)^2 + 2*P2ce*(aI^2 + aP^2)));
        
        if eta>=P1ce/P2ce+aI
            nou = eta;
        else
            nou = (gama_R*(2 * gama_2R + 1) - P1ce*(2 * gama_R + 1))...
                /(gama_R*(P2ce+2*(gama_2R+1)));
        end
        rho = max((P1ce/gama_R-1)/P2ce,nou);
        % ==============================================
        % ==============================================
        Kapa2 = 1;
        q1 = (P1ce-gama_2R+sqrt((gama_2R+1)*(P1ce^2*(1-Kapa2^2)...
            +(gama_2R-2*P1ce)*Kapa2^2)))/(G(cnt)*(-1+(gama_2R+1)*(1-Kapa2^2)));
        clear Kapa2
        
        if a12(cnt) > max(rho,l)
            if imag(q1)~=0 || isnan(q1)
                kapa1 = 1;
                p2 = P2;
                Kapa  = 0:0.001:0.999;
                Q1 = (2*P1ce/gama_2R-1)./(G(cnt)*(1-Kapa));
                kapa2 = Kapa(find(Q1>=p2, 1 ));               
                if isempty(kapa2)
                    kapa2=1;
                end
                kapa1 = min( p2 * kapa2 * G(cnt)/P1ce,1);
                if kapa1<1
                    Kapa  = 0:0.001:1;
                    Q1 = (P1ce-gama_2R+sqrt((gama_2R+1)*(P1ce^2*(1-Kapa.^2)...
                        +(gama_2R-2*P1ce)*Kapa.^2)))./(G(cnt)*(-1+(gama_2R+1)*(1-Kapa.^2)));
                    kapa2 = Kapa(find(Q1>=P2, 1 ));
                    p2=P2;
                    if isempty(kapa2)
                        kapa2=1;
                        if isreal(Q1(end))==1
                            p2=Q1(end);
                        end
                    end                    
                end                    
            elseif q1 <= P2
                kapa2 = 1;
                p2=q1;
            else
                Kapa  = 0:0.001:1;
                Q1 = (P1ce-gama_2R+sqrt((gama_2R+1)*(P1ce^2*(1-Kapa.^2)...
                    +(gama_2R-2*P1ce)*Kapa.^2)))./(G(cnt)*(-1+(gama_2R+1)*(1-Kapa.^2)));
                kapa2 = Kapa(find(Q1>=P2, 1 ));
                p2=P2;
            end
        else
            kapa2=0;
            p2=min((P1ce-gama_2R+sqrt((gama_2R+1)*(P1ce^2*(1-kapa2^2)...
            +(gama_2R-2*P1ce)*kapa2^2)))/(G(cnt)*(-1+(gama_2R+1)*(1-kapa2^2))),P2);
         end
        
%         p2 = min(q_kapa2,P2ce);
        kapa1 = min( p2 * kapa2 * G(cnt)/P1ce,1);
        
        % ==============================================
        % ==============================================
        
        phi2_a = 0; phi1_a = pi;
        % ===========================================
        % ======== Worst-case channels ==============
        % ===========================================
        phi2 = phi2_a - d_ch1;
        phi1 = phi1_a;
        clear pp2 
        if isreal(p2)==1
            pp2=p2:-0.001:0;
        else
            keyboard;
        end
        

        
        
        P1cr = P1*H(cnt)*(1-beta)^2;
        cte  = (pp2 * G(cnt)*(1+beta)^2 + 1).^2 - abs(pp2...
            *kapa2*G(cnt)*(1+beta)^2).^2;
        
        R11_imp_non = .5 * log2(((pp2*G(cnt)*(1+beta)^2 + P1cr + 1).^2)./(cte)...
            -(abs(pp2*exp(1i*phi2)*kapa2*G(cnt)*(1+beta)^2 + P1cr*exp(1i*(phi1))*kapa1).^2)./(cte));
        
        
        
        if isempty(find(R11_imp_non>=R1_imp_th(Asho2,bsho), 1 ))
            p2=0;
            R1_imp_non(Asho2,bsho) = R11_imp_non(1);
        else
            p2=pp2(find(R11_imp_non>=R1_imp_th(Asho2,bsho), 1 ));
            R1_imp_non(Asho2,bsho) = R11_imp_non(find(R11_imp_non>=R1_imp_th(Asho2,bsho), 1 ));
        end
        
        R2_impth_non(Asho2,bsho) = .5 * log2(1+p2*F(cnt)*((1-beta)^2)*...
            (F(cnt)*p2*(1-kapa2^2)*(1-beta)^2+2));
        % ===============================================
        % ============ Proper Signaling =================
        % ===============================================
        
        p2_pr   = min((P1ce/(2^Rnon-1)-1)/G(cnt),P2);
        clear pp2
        pp2=p2_pr:-0.001:0;
        %p2_pr_r = p2_pr*(abs(f_r(cnt))/abs(f(cnt)))^2;
        
        R11_p_non = log2(1+P1cr./(1+pp2*G(cnt)*(1+beta)^2));
        
        if isempty(find(R11_p_non>=R1_p_th(Asho2,bsho), 1 ))
            p2_pr=0;
            R1_p_non(Asho2,bsho) = R11_p_non(1);
        else
            p2_pr=pp2(find(R11_p_non>=R1_p_th(Asho2,bsho), 1 ));
            R1_p_non(Asho2,bsho) = R11_p_non(find(R11_p_non>=R1_p_th(Asho2,bsho), 1 ));
        end
        
        
%         R2_p_non(Asho2,bsho) = log2(1+p2_pr*abs(f_r(cnt))^2);
        R2_pth_non(Asho2,bsho) = log2(1+p2_pr*F(cnt)*(1-beta)^2);
        
     
    end
    
    ave_impth2_2IC (Asho2) = sum(R2_impth_2IC(Asho2,:))/iter;
    ave_impth2 (Asho2) = sum(R2_impth(Asho2,:))/iter;
    ave2_pth2 (Asho2)  = sum(R2_pth(Asho2,:))/iter;
    bsho = 0;
    
    ave_impth2_non (Asho2) = sum(R2_impth_non(Asho2,:))/iter;
    ave2_pth2_non (Asho2)  = sum(R2_pth_non(Asho2,:))/iter;
       
    R1_impth(Asho2)=sum(R1_imp_th(Asho2,:))/iter;
    R1_pth(Asho2)=sum(R1_p_th(Asho2,:))/iter;
    
    Rsum_impth(Asho2) = R1_impth(Asho2) + ave_impth2(Asho2);
    Rsum_pth(Asho2)   = R1_pth(Asho2)   + ave2_pth2(Asho2);
    
    R1_impth_non(Asho2) = sum(R1_imp_non(Asho2,:))/iter;
    R1_pth_non(Asho2) = sum(R1_p_non(Asho2,:))/iter;
    
    Rsum_impth_non(Asho2) =  R1_impth_non(Asho2) + ave_impth2_non(Asho2);
    Rsum_pth_non(Asho2)   = R1_pth_non(Asho2)   + ave2_pth2_non(Asho2);
    
end

eps2=10*log10(eps2);
figure('DefaultAxesFontSize',40)
hold on
plot(eps2,ave_impth2,'LineWidth',2)
plot(eps2,ave_impth2_2IC,'LineWidth',2)
%plot(0:0.01:0.1,ave2_pth2,':','LineWidth',8)
plot(eps2,ave2_pth2_non,':','LineWidth',2)
plot(eps2,ave_impth2_non,'LineWidth',2)
grid on
xlabel('SNR')
ylabel('Average Rate  of User 2')
legend('R-IGS','R-IGS-2IC','R-PGS','N-IGS')
legend('Location','best')


figure('DefaultAxesFontSize',40)
hold on
plot(eps2,100*(ave_impth2-ave_impth2_2IC)./ave_impth2_2IC,'LineWidth',2)
grid on
xlabel('SNR')





