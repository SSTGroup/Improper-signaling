% --------------------------------------------------------
% --------------------------------------------------------
% In this code, we obtain the achievable rate region by our main algorithm.
% It can simply be used for deriving other figures in the paper

clc,
clear;
close all;
% ==========================================
% =========== Initialization ===============
global P1 P2 sigma2 sigma_eta sigma_eta_c alpha
P1=10; % The maximum total power of user 1
P2=10; % The maximum total power of user 2
sigma2 = 1;
sigm=0.;
sigma_eta=0;
sigma_eta_c = sigm*sigma_eta;
% ==========================================
iter=1;
RandStream.setGlobalStream(RandStream('mcg16807','seed',sum(100*clock)));
% ==========================================
% ========== Estimated Channels ============
% ==========================================
h11 = (1/sqrt(2))*(randn(1,iter)+1i*randn(1,iter)); 

h12 = (1/sqrt(2))*(randn(1,iter)+1i*randn(1,iter));

h21 = (1/sqrt(2))*(randn(1,iter)+1i*randn(1,iter)); 

h22 = (1/sqrt(2))*(randn(1,iter)+1i*randn(1,iter)); 
% --------------------------------------------------



L=1000;
% alpha1=0:0.2:0.6;
dummy=0; %alpha=0.5;


eps2=[0,0.1,0.3,0.5];
alpha1=0:0.01:1;

% R1_FP=zeros(length(eps2),iter);
for alpha=alpha1%eps1=eps2
    
    dummy = dummy+1;
    for cnt=1:iter
       
        h=[h11(cnt),h12(cnt);h21(cnt),h22(cnt)];
        
        % -----------------------------------------------------
        % -----------------------------------------------------
        % The first algorithm in the paper, which is based on fractional
        % programming
        [R1_FP(dummy,cnt),R2_FP(dummy,cnt),PV_FP(:,dummy,cnt), ...
            P_FP(:,dummy,cnt),time_FP(dummy,cnt)] = IGS_FP(h);
        
        
        % The smplified algorithm in the paper
        [R1_i_sub(dummy,cnt),R2_i_sub(dummy,cnt),R1_p_sub(dummy,cnt), ...
            R2_p_sub(dummy,cnt),kap_1(dummy,cnt),kap_1(dummy,cnt), ...
            P_sub_1(dummy,cnt),P_sub_1(dummy,cnt),time_sub(dummy,cnt)] ...
            = IGS_Dis(h);
        
                
        % -----------------------------------------------------
        % -----------------------------------------------------
        
         
        
        save data 
       
        
    end
    R_FP_1(dummy) = mean(R1_FP(dummy,:));
    R_FP_2(dummy) = mean(R2_FP(dummy,:));
    R_sub_1(dummy) = mean(R1_i_sub(dummy,:));
    R_sub_2(dummy) = mean(R2_i_sub(dummy,:));
    t_FP(dummy)=mean(time_FP(dummy,:));
    t_sub(dummy)=mean(time_sub(dummy,:));
    R_p_1(dummy) = mean(R1_p_sub(dummy,:));
    R_p_2(dummy) = mean(R2_p_sub(dummy,:));

    
end

figure('DefaultAxesFontSize',40)
hold on
plot(R_FP_1,R_FP_2,'LineWidth',8)
plot(R_sub_1,R_sub_2,'LineWidth',8)
plot(R_p_1,R_p_2,'LineWidth',8)
grid on
xlabel('Rate of user 1')
ylabel('Rate of user 2')
legend('IGS-FP','IGS-S','PGS')