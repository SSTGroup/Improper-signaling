%% Ergodic capacity curves for simmetric channels with maximum power
% Ergodic Rate as function of inr
clear %all
% close all

SNRdB = 15; SNR = 10^(SNRdB/10);
INRdB = [0:2:6,6.3,(SNRdB-1),(SNRdB+1):2:30,30]; INR = 10.^(INRdB/10);

% channels at Rx
snr = SNR;
%inr = INR;
k=1;

r_pp=log2(exp(1))*exp(1/snr(1))*expint(1/snr(1));
%% PGS

alpha1 = 1;              % Maximum power for user 1
alpha2 = 1;              % Maximum power for user 2

RaPGS = [];
RaJPGS = [];  % Jensens app.

for nn = 1:length(INR)
    %alpha2 = alpha(nn);
    inr = INR(nn);
    if snr*alpha1 ~= inr*alpha2
        RaPGS = [RaPGS log2(exp(1))*( snr*alpha1./(snr*alpha1-inr*alpha2)).*....
            (exp(1./(alpha1*snr)).*expint(1./(alpha1*snr)) - exp(1./(inr*alpha2)).*expint(1./(alpha2*inr)))]; 
%     else
%         RaPGS = [RaPGS log2(exp(1)).*exp(1./(alpha1*snr)).* ...
%             expint(2,1./(alpha1*snr))];
    end
    RaJPGS = [RaJPGS log2(1+((snr.*alpha1)/(inr.*alpha2)).*exp(1./(inr.*alpha2))...
        .*expint(1./(inr.*alpha2)))];
end

%% IGS - 1 user

alpha1 = 1;              % Maximum power for user 1
alpha2 = 1;              % Maximum power for user 2

RaIGS1_1 = [];
RaIGS1_2 = [];  % Jensens app.

for nn = 1:length(INR)
    %alpha2 = alpha(nn);
    inr = INR(nn);
    if snr*alpha1 ~= inr*alpha2
        RaIGS1_1 = [RaIGS1_1 0.5*log2(exp(1))*( exp(1./(alpha1*snr)).* ...
            expint(1./(alpha1*snr)) + ( snr*alpha1./(snr*alpha1-2*inr*alpha2)).*....
            (exp(1./(alpha1*snr)).*expint(1./(alpha1*snr)) - ...
            exp(1./(2*inr*alpha2)).*expint(1./(2*alpha2*inr))))]; 
        
        RaIGS1_2 = [RaIGS1_2 0.5*log2(exp(1))*( 2*snr*alpha1./(2*snr*alpha1-inr*alpha2)).*....
            (exp(1./(2*alpha1*snr)).*expint(1./(2*alpha1*snr)) - ...
            exp(1./(inr*alpha2)).*expint(1./(alpha2*inr)))]; 
    end
end


%% IGS - 1 user- Exhaustive search

for nn=1:length(INR)
    [r_2IGS(nn),kap(nn)]=PGS_IGS_Exhaustive_kappa(snr,INR(nn));
end





figure; hold on
plot(INRdB,RaIGS1_1+RaIGS1_2);
plot(INRdB,2*RaPGS);
plot(INRdB,r_2IGS,'*')
grid on
title('Ergodic Capacity');
legend('IGS','PGS','IGS-E');
xlabel('INR (dB)');ylabel('Sum Rate (bps/Hz)');hold off;



save data

