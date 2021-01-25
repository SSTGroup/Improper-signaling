function [R1_FP,R2_FP,PV_FP,P_FP,time_FP]=IGS_FP(h)

global P1 P2 sigma2 sigma_eta sigma_eta_c alpha
% -----------------------------------------------
% --------------- Coefficients ------------------
% -----------------------------------------------
a1 = [abs(h(1,1))^2*(1+sigma_eta);...
    abs(h(2,1))^2*(1+sigma_eta)];
a2 = [abs(h(1,2))^2*(1+sigma_eta);...
    abs(h(2,2))^2*(1+sigma_eta)];

b1 = [abs(h(1,1))^2*(sigma_eta);...
    abs(h(2,1))^2*(1+sigma_eta)];
b2 = [abs(h(1,2))^2*(1+sigma_eta);...
    abs(h(2,2))^2*(sigma_eta)];

f1 = [h(1,1)'^2 h(2,1)'^2]'; 
f2 = [h(1,2)'^2 h(2,2)'^2]'; 

f1_til = [h(1,1)'^2*(sigma_eta_c)'...
    h(2,1)'^2*(sigma_eta_c)']'; 
f2_til = [h(1,2)'^2*(sigma_eta_c)'...
    h(2,2)'^2*(sigma_eta_c)']'; 

g1 = [0 h(2,1)'^2]'; 
g2 = [h(1,2)'^2 0]'; 

% -----------------------------------------------
% -----------------------------------------------

tic
if alpha==0 
   P_FP = [0 P2];
    PV_FP = [0 -sigma_eta]*P2;
    
    
    R1_FP = 0;
    R2_FP = rate_IGS(sigma2,P_FP,PV_FP,a2,b2,f2,g2,f2_til);
    
elseif alpha==1
    P_FP = [P1 0];
    PV_FP = [-sigma_eta 0]*P1;
    
    R1_FP = rate_IGS(sigma2,P_FP,PV_FP,a1,b1,f1,g1,f1_til);
    R2_FP = 0;
else
    p_l(1) = min(2*alpha,1)*P1;
    p_l(2) = min(2*(1-alpha),1)*P2;
    
    q_l = p_l;%[0 0];
    
    cnstrt_2=1;count_2=0;
    
    while cnstrt_2 == 1&&count_2<20
        count_2 = count_2+1;
        cnstrt  = 1;count=0;cvx_clear;%r_pre=0;
        
        %-----------------------------------------------
        temp1=((sigma2^2+p_l*a1)^2-abs(q_l*f1+p_l*f1_til)^2) ...
        /((sigma2^2+p_l*b1)^2-abs(q_l*g1+p_l*f1_til)^2);

        temp2=((sigma2^2+p_l*a2)^2-abs(q_l*f2+p_l*f2_til)^2) ...
            /((sigma2^2+p_l*b2)^2-abs(q_l*g2+p_l*f2_til)^2);

        mu = min(temp1/alpha,temp2/(1-alpha));
        mu_l=mu;
        clear temp1 temp2
        %-----------------------------------------------
        
        while cnstrt == 1 && count<40
            
            count=count+1;
            cvx_begin sdp quiet
               variable c(1,2) nonnegative
               variable q(1,2) complex
               variable r

               temp1 = (sigma2^2+p_l*a1)^2+2*(sigma2^2+p_l*a1)* ...
                   (c-p_l)*a1;

               temp2 = (sigma2^2+p_l*a2)^2+2*(sigma2^2+p_l*a2)* ...
                   (c-p_l)*a2;

               tp1=abs(q_l*g1+p_l*f1_til)^2 + 2*real((q_l*g1+ ... 
                   p_l*f1_til)'*(c-p_l)*f1_til) ...
                   +2*real((q_l*g1+p_l*f1_til)'*(q-q_l)*g1);

               tp2=abs(q_l*g2+p_l*f2_til)^2+2*real((q_l*g2+ ...
                   p_l*f2_til)'*(c-p_l)*f2_til) ...
                   +2*real((q_l*g2+p_l*f2_til)'*(q-q_l)*g2);

               maximize r
               subject to
                 c(1) <= P1;
                 c(2) <= P2;
                 abs(q(1)) <= c(1);
                 abs(q(2)) <= c(2);



                -alpha*mu*(sigma2^2+c*b1)^2 - (q*f1+c*f1_til)*(q*f1+c*f1_til)' ...
                    +temp1 + alpha*mu*tp1>=r;

                -(1-alpha)*mu*(sigma2^2+c*b2)^2 - (q*f2+c*f2_til)*(q*f2+c*f2_til)' ...
                    +temp2 + (1-alpha)*mu*tp2>=r;

            cvx_end
            
            temp3=((sigma2^2+c*a1)^2-abs(q*f1+c*f1_til)^2)- ...
                mu*alpha*((sigma2^2+c*b1)^2-abs(q*g1+c*f1_til)^2);

            temp4=((sigma2^2+c*a2)^2-abs(q*f2+c*f2_til)^2)- ...
                mu*(1-alpha)*((sigma2^2+c*b2)^2-abs(q*g2+c*f2_til)^2);

            if min(temp3,temp4)<=0.01|| mu_l==mu
                cnstrt=0;
                 mu_l=mu;
            else
                %-----------------------------------------------
                temp1=((sigma2^2+c*a1)^2-abs(q*f1+c*f1_til)^2) ...
                /((sigma2^2+c*b1)^2-abs(q*g1+c*f1_til)^2);

                temp2=((sigma2^2+c*a2)^2-abs(q*f2+c*f2_til)^2) ...
                    /((sigma2^2+c*b2)^2-abs(q*g2+c*f2_til)^2);

                mu = min(temp1/alpha,temp2/(1-alpha));
                clear temp1 temp2
                %-----------------------------------------------
                
            end          

        end
        
         r1_l=rate_IGS(sigma2,p_l,q_l,a1,b1,f1,g1,f1_til);
        r2_l=rate_IGS(sigma2,p_l,q_l,a2,b2,f2,g2,f2_til);
        
        r1_ln=rate_IGS(sigma2,c,q,a1,b1,f1,g1,f1_til);
        r2_ln=rate_IGS(sigma2,c,q,a2,b2,f2,g2,f2_til);
        
        if abs((r1_ln-r1_l))/r1_l<0.0001&&abs((r2_ln-r2_l))/r2_l<0.0001
            cnstrt_2 = 0;
          
        end 
        p_l = c;
        q_l = q;

        
        
%         if sum(abs((c-p_l))./p_l)<0.01&&sum(abs((q-q_l)./p_l)) < 0.01
%             cnstrt_2 = 0;            
%         else
%             r_pre=r;
%         end  
%         
%         p_l=c; q_l=q;


    end
    P_FP = p_l;
    PV_FP = q_l;
%     temp1=((sigma2^2+c_l*a1)^2-abs(q_l*f1+c_l*f1_til)^2) ...
%         /((sigma2^2+c_l*b1)^2-abs(q_l*g1+c_l*f1_til)^2);
%     temp2=((sigma2^2+c_l*a2)^2-abs(q_l*f2+c_l*f2_til)^2) ...
%         /((sigma2^2+c_l*b2)^2-abs(q_l*g2+c_l*f2_til)^2);
    R1_FP=rate_IGS(sigma2,p_l,q_l,a1,b1,f1,g1,f1_til);
    R2_FP=rate_IGS(sigma2,p_l,q_l,a2,b2,f2,g2,f2_til);
    
end
time_FP=toc;