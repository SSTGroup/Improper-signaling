tic
clc,
clear;
% close all;
% ==========================================
% =========== Initialization ===============
P=50; % The maximum total power of users 
sigma2 = 1; % Noise Power
K=5; % Number of Users (K)
% ==========================================
iter=100;
RandStream.setGlobalStream(RandStream('mcg16807','seed',sum(100*clock)));
% ==========================================
% ========= Generating Channels ============
% ==========================================
h = (1/sqrt(2))*(randn(K,K,iter)+1i*randn(K,K,iter));  
H = abs(h).^2; % Channel gains
% --------------------------------------------------------------
% --------------------------------------------------------------



eta=3;
gamma=zeros(1,K);
dummy=0; 
alpha=(1/K)*ones(1,K);
P_c2=40;%[10,20,40,60,80];
alpha2=min(alpha);
for P_c=P_c2
    dummy = dummy+1;
    for cnt=1:iter
        p_pre = P*ones(1,K); % Initial Powers
        p_new = p_pre;

        for k=1:K                
            temp = sigma2^2 + sum(p_pre.*H(k,:,cnt))-p_pre(k)*H(k,k,cnt);                
            gamma(k) = (p_pre(k)*H(k,k,cnt))/temp;                
        end
        % ===========================================
        % ===========================================
        u = log2(1 + gamma)./(alpha.*(eta.*p_pre+P_c));
        mu = min(u);
        % ----------------------------------------------
        % ----------------------------------------------
        cnstrt=1; iter2=0;
        while cnstrt==1 && iter2<40
            iter2 = iter2+1;
            p_pre = p_new; 
            % ===========================================
            % ===========================================
            a = gamma./(1+gamma);
            % ===========================================
            % ===========================================
            b = log2(1+gamma) - (gamma.*log2(gamma))./(1+gamma);
            % ===========================================
            % =========================================== 
            cvx_begin sdp quiet
               variable q(1,K) %nonnegative
               variable u nonnegative
               clear tp


               tp(1)=(2^q(2)) * H(1,2,cnt);
               if K>2
                   for j=3:K
                       tp(1) = (2^q(j)) * H(1,j,cnt) + tp(1);
                   end
               end
               for k=2:K
                   tp(k)=0;
                   for j=1:K
                       if j~=k
                            tp(k) = (2^q(j)) * H(k,j,cnt) + tp(k);
                       end
                   end
               end

               maximize u
               subject to  
                for k=1:K
                    q(k) <= log2(P);
                    a(k)*q(k)+a(k)*log2(H(k,k,cnt)) - ...
                        (a(k)/log(2)) * log(sigma2^2 +tp(k)) + b(k) - ...
                     mu*alpha(k)*(eta*2^q(k)+P_c)>=u;
                end   

            cvx_end
            % ============================================
            % ============================================
            p_new = 2.^q;

            for k=1:K                
                temp = sigma2^2 + sum(p_new.*H(k,:,cnt))-p_new(k)*H(k,k,cnt);                
                gamma(k) = (p_new(k)*H(k,k,cnt))/temp;                
            end

            u1 = log2(1 + gamma)- mu.*(alpha.*(eta.*p_new+P_c));
            t = min(u1);

            u = log2(1 + gamma)./(alpha.*(eta.*p_new+P_c));
            mu_new = min(u);


            % ============================================
            % ============================================
            if t>10^(-4)
                mu=mu_new;
            elseif t<10^(-4) || t>=0
                cnstrt=0;
                mu_opt=mu_new;
                p_opt=p_new;
            else%if t<0
                cnstrt=0;
                mu_opt=mu;
                p_opt=p_pre;
            end
            if iter2==40
                mu_opt = mu_new;
                p_opt  = p_new;
            end
        end
        EE_p (dummy,cnt)  = mu_opt;
        % --------------------------------------
        p1_opt(:,dummy,cnt) = p_opt;
        % --------------------------------------
        R1_p_sub(:,dummy,cnt) = log2(1 + gamma);
        % ============================================
        % ========= Complementary Variances ==========
        % ============================================

        % ==================================================
        % ===== Initialization - Optimization Problem ======
        % ==================================================
        clear a
        a=H(:,:,cnt);
        % -----------------------
        clear b
        b=H(:,:,cnt);
        for k=1:K
            b(k,k)=0;
        end
        % ---------------------
        f=h(:,:,cnt).^2;
        g=f;
        for k=1:K
           g(k,k)=0; 
        end



        % -------------------------------------------------  
        % -------------------------------------------------  
        for k=1:K
            c_l(k) = p_opt(k);
        end
        % -------------------------------------------------
        % -------------------------------------------------  
        q_l = c_l;
        t_l = EE_p (dummy,cnt);
        q_ln = q_l; imp_slv=0;
        % =========================================  
        % =========================================
        cnstrt=1; count_2=0;
        while cnstrt == 1 && count_2<60
            count_2=count_2+1;
            t_u = 5*t_l;
            % --------------------------------------
            while (t_u-t_l)/t_u>10^(-3)
                t_n =(t_u+t_l)/2;
                e = 2.^(2.*alpha.*t_n.*(eta.*p_opt+P_c));
                % --------------------------------------
                cvx_begin sdp quiet
                   variable q(1,K) complex
                   clear tp
                   for k=1:K
                       tp(k) = abs(q_l*g(k,:)')^2  ...
                           +2*real((q_l*g(k,:)')'*(q-q_l)*g(k,:)');
                   end


                   maximize 1
                   subject to
                       for k=1:K
                           abs(q(k)) <= c_l(k);
                           (sigma2^2+c_l*a(k,:)')^2 - e(k)*(sigma2^2+c_l*b(k,:)')^2 - ...
                               (q*f(k,:)')*(q*f(k,:)')' + e(k)*tp(k) >=0;
                       end

                cvx_end
                % --------------------------------------
                if strcmpi(cvx_status,'solved')
                    t_l = t_n;
                    q_ln = q; imp_slv=1;
                elseif strcmpi(cvx_status,'infeasible') %|| strcmpi(cvx_status,'failed')
                    t_u = t_n;
                else
    %                 keyboard;
                    t_l = t_n; imp_slv=1;
                    q_ln=q;
                end



            end

            % =================================
            % =================================
            if sum(abs((q_ln-q_l)./c_l)) < 0.005%0.01
                cnstrt = 0;
            end   
            q_l=q_ln;
        end

        if imp_slv==1
            % ----------------------------------------------
            for k=1:K
                clear temp
                temp=((sigma2^2+c_l*a(k,:)')^2-abs(q_l*f(k,:)')^2) ...
                /((sigma2^2+c_l*b(k,:)')^2-abs(q_l*g(k,:)')^2);
                R_i(k,dummy,cnt)=0.5*log2(temp);
                EE_i_k(k,dummy,cnt)=R_i(k,dummy,cnt)/ ...
                    (alpha(k)*(eta*p_opt(k)+P_c));

            end
            EE_i(dummy,cnt)=min(EE_i_k(:,dummy,cnt));
        else
            % --------------------------------------
            EE_i(dummy,cnt) = EE_p(dummy,cnt);
        end
                 

            
           
        
        
  
    end

    % ===========================================
    % ===========================================
    EE_i_a (dummy) = sum(EE_i(dummy,:))/iter;
    EE_p_a (dummy) = sum(EE_p(dummy,:))/iter;
    
    imp (dummy) = sum((EE_i(dummy,:)-EE_p(dummy,:))./EE_p(dummy,:))/iter;
    save data

end

