clc,
clear;
close all;
% ==========================================
% =========== Initialization ===============
global P1 P2 sigma sigma_d sigma_r alpha N_t N_r K coun coun2 coun_i coun2_i
P1=10; % The maximum total power of user 1
P2=10; % The maximum total power of user 2
sigma = 1;

sigma_d=0.2;
sigma_r=0.;
N_t = 2; N_r = 2;
K=2;

% ==========================================
iter=3;
RandStream.setGlobalStream(RandStream('mcg16807','seed',sum(100*clock)));
% ==========================================
% ========== Generating Channels ============
% ==========================================
h = (1/sqrt(2))*(randn(K*N_r,K*N_t,iter)+1i*randn(K*N_r,K*N_t,iter));  
% --------------------------------------------------
coun=0;
coun2=0;
coun_i=0;
coun2_i=0;

for cnt=1:iter
    cnt
    dum=0;
    for alpha=0:.05:1
        dum=dum+1
      
        % ----------------------------------------------------
        p_new_1 = (P1/(N_t))*eye(N_t,N_t); 
        p_new_2 = (P2/(N_t))*eye(N_t,N_t); 
        % ----------------------------------------------------      
        [r_uni(cnt,dum), r_u(cnt,dum), R_1_u(cnt,dum), R_2_u(cnt,dum), p_new_1, ...
            p_new_2, t_real_u(cnt)]= r_MIMO_1(h(:,:,cnt), p_new_1, p_new_2);
        
%         eval(sprintf('save rr_Nt_%i_Nr_%i_sigma_%i',N_t,N_r,10*sigma_d));
        
    end
    figure; plot(sort(R_1_u(cnt,:),'descend'),sort(R_2_u(cnt,:)))
    hold on
end
save rr_2

