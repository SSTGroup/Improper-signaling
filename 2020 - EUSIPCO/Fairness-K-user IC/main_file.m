clc,
clear;
% ==========================================
% =========== Initialization ===============
global P P1 P2 sigma sigma_d sigma_r alpha N_t N_r K coun coun2 coun_i coun2_i
sigma = 1;
P=10;
% sigma_d=0.5;
N_t = 6; N_r = 6;
K=3;
alpha=1/K;
% ==========================================
iter=100;
RandStream.setGlobalStream(RandStream('mcg16807','seed',sum(100*clock)));
% % ==========================================
% % ========== Estimated Channels ============
% % ==========================================
h = (1/sqrt(2))*(randn(K*N_r,K*N_t,iter)+1i*randn(K*N_r,K*N_t,iter));  % Secondary Channel Gain
% % --------------------------------------------------
coun=0;
coun2=0;
coun_i=0;
coun2_i=0;

e1=[0,0.2,.4,.6,.8,1];
for cnt=1:iter
    for dum=1:length(e1)
        sigma_d=e1(dum);

        [~, r(cnt,dum), ~, t_(cnt,dum), P_n] ...
            = r_MIMO_k(h(:,:,cnt));
        
        
        [r_i(cnt,dum), t_i(cnt,dum)] ...
            = r_MIMO_i_k(h(:,:,cnt));
        
        
        
%         % ==============================
%         % Extracting Fig 6
%         if dum==1
%             p_new=P_n;
%             r_i(cnt,dum)=r(cnt,dum);
%         else
%             %[r_i(cnt,dum), t_i(cnt,dum)] ...
%             r_i(cnt,dum) = rate_comp(h(:,:,cnt),p_new);
%         end
%         % ===========================

        
        
        (r-r_i)./r_i
        (mean(r)-mean(r_i))./mean(r_i)
        % ------------------------------------------
        eval(sprintf('save r_Nt_%i_Nr_%i_K_%i_P_%i_vs_sigma_2',N_t,N_r,K,P));  

    end
end
