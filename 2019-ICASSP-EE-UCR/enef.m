
clear all
% close all
iter=1;iter2=1;%10000;
sigma2 =1;
kapa  =0:.01:1;
% P=[1, 5, 10, 20, 50, 70, 100];%
p  =20;
% P_S=inf;
P_S=20;
alpha =0.7; 
zeta=1/0.35;
% P_C=0.7*P;%
p_c = 2;
p_c=p_c/zeta;
h11=(1/sqrt(2))*(randn(1,iter2)+1i*randn(1,iter2));
h12=(1/sqrt(2))*(randn(1,iter2)+1i*randn(1,iter2));
h21=(1/sqrt(2))*(randn(1,iter2)+1i*randn(1,iter2));
h22=(1/sqrt(2))*(randn(1,iter2)+1i*randn(1,iter2));
dummy=0;eps=0.1:.01:2;
for eps2=eps%p=P
    dummy=dummy+1
%     p_c=P_C(dummy);
    for cnt=1:iter2
        d=1;%h12(cnt);
        f=1;%h22(cnt);
        g=eps2;%h21(cnt);
        h=1;%h11(cnt);


        Rmax  = log2(1+p*abs(h)^2/sigma2);
        gamma_1  = p*abs(h)^2/sigma2;
        gamma_2  = 2^(2*Rmax)-1;
        gamma_2a = 2^(2*alpha*Rmax)-1;
        clear U q
        beta = 1-gamma_1/gamma_2a;
        if beta<(abs(g)^2*(sigma2+p*abs(d)^2))/(sigma2*abs(f)^2) 
            temp3 = sqrt(beta^2+(1-kapa.^2)*(gamma_2/gamma_2a-1))-beta;
            temp4=sigma2./(abs(g)^2*(1-kapa.^2));

            q=temp3.*temp4;

            q_p = 0:.01:q(1);
            R_p = log2(1+q_p*abs(f)^2/((sigma2+p*abs(d)^2)));
            u_p=R_p./(q_p+p_c);


            temp1 = (2*abs(f)^2/(sigma2+p*abs(d)^2))*(1-beta*(abs(f)^2)* ...
                sigma2/(abs(g)^2*(sigma2+p*abs(d)^2)));
            temp2 = (abs(f)^4*sigma2^2/(abs(g)^4*(sigma2+p*abs(d)^2)^2))* ...
                (gamma_2/gamma_2a-1);
            R = 0.5*log2(temp1*q + temp2 + 1); 

            U = R./(q+p_c);
%             U_p(dummy,cnt) = max(U(1),max(u_p));
%             U_imp(dummy,cnt) = max(U_p(dummy,cnt), max(U));
%             benef(dummy,cnt) = (max(U)-U(1))/U(1);
%             benef2(dummy,cnt)= (U_imp(dummy,cnt)-U_p(dummy,cnt))/U_p(dummy,cnt);
%             %[(U(end-1)-U(1))/U(1), (max(U)-U(1))/U(1)]
%         %     r=log2(1+q(1)*abs(f)^2/(1+p*abs(d)^2));
%     %         temp1*(q(1)+p_c)-2*log(2)*R(1)*(2^(2*R(1)))>0%(temp1(1)*q(1)+ temp2 + 1)%>-0.01
%             if temp1*(q(1)+p_c)-2*log(2)*R(1)*(2^(2*R(1)))<0 && benef(dummy,cnt)>0.001
%                 keyboard;
%             end


            % ==========================================
            % =============== P_pgs ====================
            % ==========================================


            % ========== Optimal Power-Proper =========

    %         p_l=q(1);
            mu=U(1); cntr=1;
            clear temp
            temp = abs(f)^2/(sigma2+p*abs(d)^2);
            while cntr==1
                p_l=(temp/(log(2)*mu)-1)/temp;
                temp_n = log2(1+temp*p_l) - mu*(p_l+p_c);
                if temp_n<=0.000001
                    cntr=0;
                else
                    mu=log2(1+temp*p_l) / (p_l+p_c);
                end
            end
            % ========================================
            if P_S>=max(q(end),q(end-1))
                if p_l<=q(1)
                    p_opt(dummy,cnt) = p_l;
                    U_p2 (dummy,cnt)   = log2(1+temp*p_l) / (p_l+p_c);
                    U_imp2 (dummy,cnt) = U_p2(dummy,cnt);
                else
                    U_p2 (dummy,cnt)   = U(1);
                    if temp1*(q(1)+p_c)-2*log(2)*R(1)*(2^(2*R(1)))>0
                        p_l2 = q(1);%max(q(end),q(end-1));
                        mu_2=0.5*log2(temp1*p_l2 + temp2 + 1) / (p_l2+p_c);%max(U(end),U(end-1));
                        cntr=1;
                        while cntr==1
                            p_l2=(temp1/(2*log(2)*mu_2)-1-temp2)/temp1;
                            temp_2 = 0.5*log2(temp1*p_l2 + temp2 + 1) - mu_2*(p_l2+p_c);
                            if temp_2<=0.000001
                                cntr=0;
                            else
                                mu_2=0.5*log2(temp1*p_l2 + temp2 + 1) / (p_l2+p_c);
                            end
                        end

                        if p_l2<=max(q(end),q(end-1))
                            R_theo = 0.5*log2(temp1*p_l2 + temp2 + 1);
                            U_imp2 (dummy,cnt) = R_theo/(p_l2+p_c);
                            p_opt(dummy,cnt) = p_l2;
                        else
                            U_imp2 (dummy,cnt) =max(U(end),U(end-1));
                            p_opt(dummy,cnt) =max(q(end),q(end-1));
                        end

                    else
                        U_p2 (dummy,cnt)   = U(1);
                        U_imp2 (dummy,cnt) = U_p2(dummy,cnt);  
                        p_opt(dummy,cnt) = q(1);
                    end


                end
            elseif P_S<q(1)
                if p_l<=P_S
                    p_opt(dummy,cnt) = p_l;
                    U_p2 (dummy,cnt)   = log2(1+temp*p_l) / (p_l+p_c);
                    U_imp2 (dummy,cnt) = U_p2(dummy,cnt);
                else
                    p_opt(dummy,cnt) = P_S;
                    U_p2 (dummy,cnt)   = log2(1+temp*P_S) / (P_S+p_c);
                    U_imp2 (dummy,cnt) = U_p2(dummy,cnt);
                end
            else
                
                if p_l<=q(1)
                    p_opt(dummy,cnt) = p_l;
                    U_p2 (dummy,cnt)   = log2(1+temp*p_l) / (p_l+p_c);
                    U_imp2 (dummy,cnt) = U_p2(dummy,cnt);
                else
                    U_p2 (dummy,cnt)   = U(1);
                    if temp1*(q(1)+p_c)-2*log(2)*R(1)*(2^(2*R(1)))>0
                        p_l2 = q(1);%max(q(end),q(end-1));
                        mu_2=0.5*log2(temp1*p_l2 + temp2 + 1) / (p_l2+p_c);%max(U(end),U(end-1));
                        cntr=1;
                        while cntr==1
                            p_l2=(temp1/(2*log(2)*mu_2)-1-temp2)/temp1;
                            temp_2 = 0.5*log2(temp1*p_l2 + temp2 + 1) - mu_2*(p_l2+p_c);
                            if temp_2<=0.000001
                                cntr=0;
                            else
                                mu_2=0.5*log2(temp1*p_l2 + temp2 + 1) / (p_l2+p_c);
                            end
                        end

                        if p_l2<=P_S%max(q(end),q(end-1))
                            R_theo = 0.5*log2(temp1*p_l2 + temp2 + 1);
                            U_imp2 (dummy,cnt) = R_theo/(p_l2+p_c);
                            p_opt(dummy,cnt) = p_l2;
                        else
                            R_theo = 0.5*log2(temp1*P_S + temp2 + 1);
                            U_imp2 (dummy,cnt) = R_theo/(P_S+p_c);
                            p_opt(dummy,cnt) = P_S;
                        end

                    else
                        U_p2 (dummy,cnt)   = U(1);
                        U_imp2 (dummy,cnt) = U_p2(dummy,cnt);  
                        p_opt(dummy,cnt) = q(1);
                    end
                end

            end
            if p_opt(dummy,cnt)>max(q)
                keyboard;
            end
            benef3(dummy,cnt)= (U_imp2(dummy,cnt)-U_p2(dummy,cnt))/U_p2(dummy,cnt);
    % %         temp1(1)-R(1)*(temp1(1)+ temp2 + 1)>-0.01
    %     %     temp1(1)-2*log(2)*r*(2^(2*r))>-0.01
    %     %     (temp1-2*log(2)*R.*2.^(2.*R))./(q.^2)
    %         % ==================================
    % %         figure
    % %         hold on
    % %         plot(kapa,U)
    % %         xlabel('\kappa')
    % %         ylabel('Energy Efficiency')
    % %         grid on
    % %         % ==================================
    % %         figure
    % %         plot(q_p,u_p)
    % %         xlabel('Power')
    % %         ylabel('Energy Efficiency Proper')
    % %         grid on
    %         % ==================================
    %         figure
    %         hold on
    %         plot([q_p,q],[u_p,U])
    %         plot([q(1),q(1)],[0,U_imp(cnt)*1.05])
    %         xlabel('Power')
    %         ylabel('Energy Efficiency')
    %         grid on
        % %     % ==================================
        %     figure
        %     plot(kapa,R)
        %     xlabel('\kappa')
        %     ylabel('Rate')
        %     grid on
        % %     % ==================================
        %     figure
        %     plot(kapa,q)
        %     xlabel('\kappa')
        %     ylabel('Q')
        %     grid on    

        else
            benef(dummy,cnt)=0;
            benef2(dummy,cnt)=0;
            benef3(dummy,cnt)=0;
            % ===============================================
            % ===============================================
            % ========== Optimal Power-Proper =========
            p_max=(p*abs(h)^2/(2^(alpha*Rmax)-1)-1)/(abs(g)^2);
            clear temp
            temp = abs(f)^2/(sigma2+p*abs(d)^2);
            if p_max==inf
                mu=0;
            else
                mu=log2(1+temp*p_max) / (p_max+p_c); 
            end
            cntr=1;
            while cntr==1
                p_l=(temp/(log(2)*mu)-1)/temp;
                temp_n = log2(1+temp*p_l) - mu*(p_l+p_c);
                if temp_n<=0.000001
                    cntr=0;
                else
                    mu=log2(1+temp*p_l) / (p_l+p_c);
                end
            end
            % ========================================
            % ========================================
            clear q
            q=min(p_max,P_S);
            if p_l<=q
                p_opt(dummy,cnt) = p_l;
                U_p2 (dummy,cnt)   = log2(1+temp*p_l) / (p_l+p_c);
                U_imp2 (dummy,cnt) = U_p2(dummy,cnt);
            else
                U_p2 (dummy,cnt)   = log2(1+temp*q) / (q+p_c);
                p_opt(dummy,cnt) = q(1);
                U_imp2 (dummy,cnt) = U_p2(dummy,cnt);
            end
            
            
        end
    end
%     benef_ave(dummy)  = 100*sum(benef(dummy,:))/iter2;
%     benef2_ave(dummy) = 100*sum(benef2(dummy,:))/iter2;
    benef3_ave(dummy) = 100*sum(benef3(dummy,:))/iter2;
    U_imp2_ave(dummy) = sum(U_imp2(dummy,:))/iter2;
    U_p2_ave(dummy) = sum(U_p2(dummy,:))/iter2;
    
end
% figure
% hold on
% plot(log10(P),U_imp2_ave)
% plot(log10(P),U_p2_ave)
% xlabel('SNR')
% ylabel('Energy Efficiency')
% grid on  
% % ==========================================================
% figure
% hold on
% % plot(log10(P),benef2_ave)
% plot(log10(P),benef3_ave)
% xlabel('SNR')
% ylabel('Improvements Percentage')
% grid on   
% ==========================================================
% figure
% hold on
% plot(eps,U_imp2_ave/zeta)
% plot(eps,U_p2_ave/zeta)
% xlabel('Gain of Interference Link')
% ylabel('Energy Efficiency')
% grid on  
% ==========================================================
% figure
hold on
% plot(log10(P),benef2_ave)
plot(eps,benef3_ave)
xlabel('Gain of Interference Link')
ylabel('Improvements Percentage')
grid on   
%[sum(U_p)' sum(U_p2)' sum(U_imp)' sum(U_imp2)']/iter2
% 100*((sum(U_imp2)-sum(U_p2))./sum(U_p2))