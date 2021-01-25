function [R1_i,R2_i,R1_p,R2_p,kappa_1,kappa_2, ...
    p_1,p_2,time_sub]=IGS_Dis(h)

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

sigma2_tx_c=sigma_eta_c;
sigma2_rx_c=0;

tic
if alpha==0
    
    % --------------- IGS -----------------
    c_l = [0 P2];
    q_l = [0 -sigma2_rx_c-sigma2_tx_c]*P2;
    temp2=((sigma2^2+c_l*a2)^2-abs(q_l*f2+c_l*f2_til)^2) ...
            /((sigma2^2+c_l*b2)^2-abs(q_l*g2+c_l*f2_til)^2);
    R1_i = 0;
    R2_i = 0.5*log2(temp2);
    
    kappa_1 = q_l(1);
    kappa_2 = q_l(2);
    
    p_1=c_l(1);
    p_2=c_l(2);
    
    % --------------- PGS -----------------
    c_l = [0 P2];
    q_l = [0 0];
    temp2=((sigma2^2+c_l*a2)^2-abs(q_l*f2+c_l*f2_til)^2) ...
            /((sigma2^2+c_l*b2)^2-abs(q_l*g2+c_l*f2_til)^2);
    R1_p = 0;
    R2_p = 0.5*log2(temp2);
    
elseif alpha==1
    
    % --------------- IGS -----------------
    c_l = [P1 0];
    q_l = [-sigma2_rx_c-sigma2_tx_c 0]*P1;
    temp1=((sigma2^2+c_l*a1)^2-abs(q_l*f1+c_l*f1_til)^2) ...
            /((sigma2^2+c_l*b1)^2-abs(q_l*g1+c_l*f1_til)^2);
    R1_i = 0.5*log2(temp1);
    R2_i = 0;
    
    kappa_1 = q_l(1);
    kappa_2 = q_l(2);
    p_1=c_l(1);
    p_2=c_l(2);
    % --------------- PGS -----------------
    c_l = [P1 0];
    q_l = [0 0];
    temp1=((sigma2^2+c_l*a1)^2-abs(q_l*f1+c_l*f1_til)^2) ...
            /((sigma2^2+c_l*b1)^2-abs(q_l*g1+c_l*f1_til)^2);
    R1_p = 0.5*log2(temp1);
    R2_p = 0;
    
else

    % --------- Optimizing Powers ---------
    Elb=0;
    c_l = [0 P2];
    q_l = [0 0]*P2;
    temp2=((sigma2^2+c_l*a2)^2-abs(q_l*f2+c_l*f2_til)^2) ...
            /((sigma2^2+c_l*b2)^2-abs(q_l*g2+c_l*f2_til)^2);
    c_l = [P1 0];
    q_l = [0 0]*P1;
    temp1=((sigma2^2+c_l*a1)^2-abs(q_l*f1+c_l*f1_til)^2) ...
            /((sigma2^2+c_l*b1)^2-abs(q_l*g1+c_l*f1_til)^2);
    Eub=sqrt(temp1^2+temp2^2);
    
    while (Eub-Elb)/Elb>=5*1e-4
        E1 = 1/2*(Eub+Elb);
        a_prime1 = a1-(alpha*E1+1).*b1;
        a_prime2 = a2-((1-alpha)*E1+1)*b2;
        a_prime = [a_prime1';a_prime2']; clear temp
        c_prime = [alpha*sigma2^2*E1 (1-alpha)*sigma2^2*E1]';

        if det(a_prime)>0 && a_prime(1,1)>0 && a_prime(2,2)>0
            p_prime=a_prime^(-1)*c_prime;
            if p_prime(1)<=P1 && p_prime(2)<=P2
                Elb=E1;
                p_1=p_prime(1);
                p_2=p_prime(2);
            else
                Eub=E1;
            end
        else
            Eub=E1;
        end
    end


    % -------------- Computing PGS rates -------------
    c_l(1) = p_1;
    c_l(2) = p_2;
    q_l = [0 0];%c_l;
    temp1=((sigma2^2+c_l*a1)^2-abs(q_l*f1+c_l*f1_til)^2) ...
        /((sigma2^2+c_l*b1)^2-abs(q_l*g1+c_l*f1_til)^2);
    temp2=((sigma2^2+c_l*a2)^2-abs(q_l*f2+c_l*f2_til)^2) ...
        /((sigma2^2+c_l*b2)^2-abs(q_l*g2+c_l*f2_til)^2);
    R1_p = 0.5*log2(temp1);
    R2_p = 0.5*log2(temp2);
    e1=temp1; e2=temp2; q_l = c_l;
    % -----------------------------------------
    % -----------------------------------------
    
    
    % ----------------------------------------------------
    % --------- Optimizing complementary variances -------  
    % ----------------------------------------------------
    cnstrt=1;count_2=0;
    while cnstrt == 1 && count_2<40
        count_2=count_2+1;
        % =================================
        % =================================
        cvx_begin sdp quiet
           variable q(1,2) complex
           variable r1
           variable r2

           tp1=abs(q_l*g1+c_l*f1_til)^2  ...
               +2*real((q_l*g1+c_l*f1_til)'*(q-q_l)*g1);

           tp2=abs(q_l*g2+c_l*f2_til)^2 ...
               +2*real((q_l*g2+c_l*f2_til)'*(q-q_l)*g2);

           t = min(r1,r2);

           maximize t
           subject to
             abs(q(1)) <= c_l(1);
             abs(q(2)) <= c_l(2);

            (sigma2^2+c_l*a1)^2 - e1*(sigma2^2+c_l*b1)^2 - ...
                (q*f1+c_l*f1_til)*(q*f1+c_l*f1_til)' + e1*tp1>=r1;

            (sigma2^2+c_l*a2)^2 - e2*(sigma2^2+c_l*b2)^2 - ...
                (q*f2+c_l*f2_til)*(q*f2+c_l*f2_til)' + e2*tp2>=r2;
        cvx_end
        % =================================
        % =================================
        if sum(abs((q-q_l)./c_l)) < 0.005%0.01
            cnstrt = 0;
            q_l=q;
        else
            q_l=q;
        end            
    end
    
    % --------------- Deriving the rates -------------
    
    if t>0
        kappa_1 = q_l(1);
        kappa_2 = q_l(2);
        temp1=((sigma2^2+c_l*a1)^2-abs(q_l*f1+c_l*f1_til)^2) ...
            /((sigma2^2+c_l*b1)^2-abs(q_l*g1+c_l*f1_til)^2);
        temp2=((sigma2^2+c_l*a2)^2-abs(q_l*f2+c_l*f2_til)^2) ...
            /((sigma2^2+c_l*b2)^2-abs(q_l*g2+c_l*f2_til)^2);
        R1_i = 0.5*log2(temp1);
        R2_i = 0.5*log2(temp2);
    else
        kappa_1 = 0;
        kappa_2 = 0;
        R1_i = R1_p;
        R2_i = R2_p;
    end
end
time_sub =toc;
