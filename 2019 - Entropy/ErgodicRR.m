%% Pareto curves for E[R1] and E[R2] with PGS and IGS
% Ergodic Rate region
clear all
close all

SNRdB = 10; SNR = 10^(SNRdB/10);
INRdB = 7; INR = 10^(INRdB/10);
% channels at Rx1
snr11 = SNR;
inr21 = INR;

% channels at Rx2
snr22 = SNR;
inr12 = INR;
alpha = 0:0.01:1;

% k=0.6;
%% Maximum power for user 1 PGS

alpha1 = 1;              % Maximum power for user 1

R1aPGS = [];
R2aPGS = [];
R1aJPGS = [];  % Jensens app.
R2aJPGS = [];

for nn = 1:length(alpha)
    alpha2 = alpha(nn);
    if alpha2 == 0
        R1aPGS = [R1aPGS log2(exp(1))*exp(1/snr11)*expint(1/snr11)];
        R2aPGS = [R2aPGS 0];
        R1aJPGS = [R1aJPGS log2(1+snr11)];
        R2aJPGS = [R2aJPGS 0];
    elseif (snr11*alpha1 ~= inr21*alpha2)&&(snr22*alpha2 ~= inr12*alpha1)
        
        R1aPGS = [R1aPGS log2(exp(1))*( snr11*alpha1./(snr11*alpha1-inr21*alpha2)).*....
            (exp(1./(alpha1*snr11)).*expint(1./(alpha1*snr11)) - exp(1./(inr21*alpha2)).*expint(1./(alpha2*inr21)))];
        
        R2aPGS = [R2aPGS log2(exp(1))*( snr22*alpha2./(snr22*alpha2-inr12*alpha1)).*....
            (exp(1./(alpha2*snr22)).*expint(1./(alpha2*snr22)) - exp(1./(inr12*alpha1)).*expint(1./(alpha1*inr12)))];
        
        R1aJPGS = [R1aJPGS log2(1+((snr11.*alpha1)/(inr21.*alpha2)).*exp(1./(inr21.*alpha2))...
            .*expint(1./(inr21.*alpha2)))];
        R2aJPGS = [R2aJPGS log2(1+((snr22.*alpha2)/(inr12.*alpha1)).*exp(1./(inr12.*alpha1))...
            .*expint(1./(inr12.*alpha1)))];
    end
end
%% Maximum power for user 2
alpha2 = 1;              

R1bPGS = [];
R2bPGS = [];

R1bJPGS = [];  % Jensen app
R2bJPGS = [];  % Jensen app
for nn = 1:length(alpha)
    alpha1 = alpha(nn);
    if alpha1 == 0
        R1bPGS = [R1bPGS 0];
        R2bPGS = [R2bPGS log2(exp(1))*exp(1/snr22)*expint(1/snr22)];
        R1bJPGS = [R1bJPGS 0];
        R2bJPGS = [R2bJPGS log2(1+snr22)];
        
    elseif (snr11*alpha1 ~= inr21*alpha2)&&(snr22*alpha2 ~= inr12*alpha1)
        
        R1bPGS = [R1bPGS log2(exp(1))*( snr11*alpha1./(snr11*alpha1-inr21*alpha2)).*....
            (exp(1./(alpha1*snr11)).*expint(1./(alpha1*snr11)) - exp(1./(inr21*alpha2)).*expint(1./(alpha2*inr21)))];
        
        R2bPGS = [R2bPGS log2(exp(1))*( snr22*alpha2./(snr22*alpha2-inr12*alpha1)).*....
            (exp(1./(alpha2*snr22)).*expint(1./(alpha2*snr22)) - exp(1./(inr12*alpha1)).*expint(1./(alpha1*inr12)))];
        
        R1bJPGS = [R1bJPGS log2(1+((snr11.*alpha1)/(inr21.*alpha2)).*exp(1./(inr21.*alpha2))...
            .*expint(1./(inr21.*alpha2)))];
        R2bJPGS = [R2bJPGS log2(1+((snr22.*alpha2)/(inr12.*alpha1)).*exp(1./(inr12.*alpha1))...
            .*expint(1./(inr12.*alpha1)))];
    end
end

R1maxPGS = R1bPGS(end);
R2maxPGS = R2bPGS(end);
MPGS = max([R1aPGS R1bPGS R2aPGS R2bPGS]);


%% IGS - 1 user - MIGS

% User 1 --> Maximum Power

alpha1 = 1;              % Maximum power for user 1


R1aIGS1  = [];
R2aIGS1  = [];  

for nn = 1:length(alpha)
    alpha2 = alpha(nn);
    if alpha2==0
        R1aIGS1 = [R1aIGS1 log2(exp(1))*exp(1/snr11)*expint(1/snr11)];
        R2aIGS1 = [R2aIGS1 0];
    elseif (snr11*alpha1 ~= 2*inr21*alpha2)&&(2*snr22*alpha2~=inr12*alpha1)
        R1aIGS1 = [R1aIGS1 0.5*log2(exp(1))*( exp(1./(alpha1*snr11)).* ...
            expint(1./(alpha1*snr11)) + ( snr11*alpha1./(snr11*alpha1-2*inr21*alpha2)).*....
            (exp(1./(alpha1*snr11)).*expint(1./(alpha1*snr11)) - ...
            exp(1./(2*inr21*alpha2)).*expint(1./(2*alpha2*inr21))))]; 
        
        
        
        R2aIGS1 = [R2aIGS1 0.5*log2(exp(1))*( 2*snr22*alpha2./(2*snr22*alpha2-inr12*alpha1)).*....
            (exp(1./(2*alpha2*snr22)).*expint(1./(2*alpha2*snr22)) - ...
            exp(1./(inr12*alpha1)).*expint(1./(alpha1*inr12)))]; 
    end
end


%% IGS - 1 user - MIGS

% User 2 --> Maximum Power

alpha2 = 1;              % Maximum power for user 2

R1bIGS1  = [];
R2bIGS1  = [];  
for nn = 1:length(alpha)
    alpha1 = alpha(nn);
    if alpha1==0
        R2bIGS1 = [R2bIGS1 0.5*log2(exp(1))*exp(1/(2*snr22))*expint(1/(2*snr22))];
        R1bIGS1 = [R1bIGS1 0];
    elseif (snr11*alpha1 ~= 2*inr21*alpha2)&&(2*snr22*alpha2~=inr12*alpha1)
        R1bIGS1 = [R1bIGS1 0.5*log2(exp(1))*( exp(1./(alpha1*snr11)).* ...
            expint(1./(alpha1*snr11)) + ( snr11*alpha1./(snr11*alpha1-2*inr21*alpha2)).*....
            (exp(1./(alpha1*snr11)).*expint(1./(alpha1*snr11)) - ...
            exp(1./(2*inr21*alpha2)).*expint(1./(2*alpha2*inr21))))]; 
        
        R2bIGS1 = [R2bIGS1 0.5*log2(exp(1))*( 2*snr22*alpha2./(2*snr22*alpha2-inr12*alpha1)).*....
            (exp(1./(2*alpha2*snr22)).*expint(1./(2*alpha2*snr22)) - ...
            exp(1./(inr12*alpha1)).*expint(1./(alpha1*inr12)))]; 
    end
end


% %% IGS - 1 user - kappa
% 
% % User 1 --> Maximum Power
% 
% alpha1 = 1;              % Maximum power for user 1
% 
% 
% R1aIGS_k  = [];
% R2aIGS_k  = [];  
% 
% for nn = 1:length(alpha)
%     alpha2 = alpha(nn);
%     if alpha2==0
%         R1aIGS_k = [R1aIGS_k log2(exp(1))*exp(1/snr11)*expint(1/snr11)];
% %         r1a(nn)=log2(exp(1))*exp(1/snr11)*expint(1/snr11);
%         R2aIGS_k = [R2aIGS_k 0];
% %         r2a(nn)=0;
%     elseif (snr11*alpha1 ~= (1+k)*inr21*alpha2)&&((1+k)*snr22*alpha2~=inr12*alpha1) ...
%             &&(snr11*alpha1 ~= (1-k)*inr21*alpha2)&&((1-k)*snr22*alpha2~=inr12*alpha1)
%         R1aIGS_k = [R1aIGS_k 0.5*log2(exp(1))*( ...
%             ( snr11*alpha1./(snr11*alpha1-(1-k)*inr21*alpha2)).*....
%             (exp(1./(alpha1*snr11)).*expint(1./(alpha1*snr11)) - ...
%             exp(1./((1-k)*inr21*alpha2)).*expint(1./((1-k)*alpha2*inr21)))+ ...
%             ( snr11*alpha1./(snr11*alpha1-(1+k)*inr21*alpha2)).*....
%             (exp(1./(alpha1*snr11)).*expint(1./(alpha1*snr11)) - ...
%             exp(1./((1+k)*inr21*alpha2)).*expint(1./((1+k)*alpha2*inr21))))]; 
%         
%         
%         
%         R2aIGS_k = [R2aIGS_k 0.5*log2(exp(1))*( ...
%             ((1+k)*snr22*alpha2./((1+k)*snr22*alpha2-inr12*alpha1)).*....
%             (exp(1./((1+k)*alpha2*snr22)).*expint(1./((1+k)*alpha2*snr22))- ...
%             exp(1./(inr12*alpha1)).*expint(1./(alpha1*inr12)))+...
%             (1-k)*snr22*alpha2./((1-k)*snr22*alpha2-inr12*alpha1).*....
%             (exp(1./((1-k)*alpha2*snr22)).*expint(1./((1-k)*alpha2*snr22)) - ...
%             exp(1./(inr12*alpha1)).*expint(1./(alpha1*inr12))))]; 
%     end
% %     r1a(nn)=R1aIGS1(nn);
% %     r2a(nn)=R2aIGS1(nn);
% end
% 
% 
% %% IGS - 1 user
% 
% % User 2 --> Maximum Power
% 
% alpha2 = 1;              % Maximum power for user 2
% % k=0.5;
% R1bIGS_k  = [];
% R2bIGS_k  = [];  
% for nn = 1:length(alpha)
%     alpha1 = alpha(nn);
%     if alpha1==0
%         R2bIGS_k = [R2bIGS_k 0.5*log2(exp(1))*(...
%             exp(1/((1+k)*snr22))*expint(1/((1+k)*snr22)) + ...
%             exp(1/((1-k)*snr22))*expint(1/((1-k)*snr22)))];
%         R1bIGS_k = [R1bIGS_k 0];
%     elseif (snr11*alpha1 ~= (1+k)*inr21*alpha2)&&((1+k)*snr22*alpha2~=inr12*alpha1) ...
%             &&(snr11*alpha1 ~= (1-k)*inr21*alpha2)&&((1-k)*snr22*alpha2~=inr12*alpha1)
%         R1bIGS_k = [R1bIGS_k 0.5*log2(exp(1))*( ...
%             ( snr11*alpha1./(snr11*alpha1-(1-k)*inr21*alpha2)).*....
%             (exp(1./(alpha1*snr11)).*expint(1./(alpha1*snr11)) - ...
%             exp(1./((1-k)*inr21*alpha2)).*expint(1./((1-k)*alpha2*inr21)))+ ...
%             ( snr11*alpha1./(snr11*alpha1-(1+k)*inr21*alpha2)).*....
%             (exp(1./(alpha1*snr11)).*expint(1./(alpha1*snr11)) - ...
%             exp(1./((1+k)*inr21*alpha2)).*expint(1./((1+k)*alpha2*inr21))))]; 
%         
%         R2bIGS_k = [R2bIGS_k 0.5*log2(exp(1))*( ...
%             ((1+k)*snr22*alpha2./((1+k)*snr22*alpha2-inr12*alpha1)).*....
%             (exp(1./((1+k)*alpha2*snr22)).*expint(1./((1+k)*alpha2*snr22))- ...
%             exp(1./(inr12*alpha1)).*expint(1./(alpha1*inr12)))+...
%             (1-k)*snr22*alpha2./((1-k)*snr22*alpha2-inr12*alpha1).*....
%             (exp(1./((1-k)*alpha2*snr22)).*expint(1./((1-k)*alpha2*snr22)) - ...
%             exp(1./(inr12*alpha1)).*expint(1./(alpha1*inr12))))];  
%     end
% %     r1b(nn)=R1bIGS1(nn);
% %     r2b(nn)=R2bIGS1(nn);
% end




%% IGS - Exhaustive Search

alpha2 = 1;              % Maximum power for user 2
% k=0.5;
R1bIGS_ran  = [];
R2bIGS_ran  = [];
N=800;
a1=rand(1,N);
a2=rand(1,N);
a3=rand(1,N);
R1bIGS_ran=zeros(1,N);
R2bIGS_ran=zeros(1,N);
for nn = 1:N
    [R1bIGS_ran(2*10^4+nn),R2bIGS_ran(2*10^4+nn)] = two_users_IGS(a1(nn)*snr11, ...
        snr22,inr12,a1(nn)*inr21,a2(nn),a3(nn),1);
end





figure;clf;hold on; 
plot([R1bIGS_ran R2bIGS_ran],[R2bIGS_ran R1bIGS_ran],'g.');
plot([R1bIGS1 fliplr(R1aIGS1)],[R2bIGS1 fliplr(R2aIGS1)],'b','LineWidth',3);
plot([R2bIGS1 fliplr(R2aIGS1)],[R1bIGS1 fliplr(R1aIGS1)],'k','LineWidth',3);
plot([R1bPGS fliplr(R1aPGS)],[R2bPGS fliplr(R2aPGS)],'r--','LineWidth',2);
legend('E-IGS', 'MIGS-U2','MIGS-U1','PGS' );
xlabel('R_1 (bps/Hz)');ylabel('R_2 (bps/Hz)');%hold off;
grid on
plot([R2aIGS1(end) R1aIGS1(end) ],[R1aIGS1(end) R2aIGS1(end)])
plot(R1maxPGS,R2maxPGS,'rs')
