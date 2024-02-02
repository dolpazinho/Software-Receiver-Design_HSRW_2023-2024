% Mystery A Receiver

clear all;
close all;
clc;

load('mysteryA.mat')
r_A = r.';


%==========> Signal parameters
M = 4; 
refC = pammod(0:M-1,M); % 4-PAM Reference Constellation
Fs_A = 700e3; %[Hz] Sampling freq.
Fc_A = 1.6e6; %[Hz] Carrier freq.
beta_A = 0.24; % rolloff parameter for srrc
Symbol_Period_A = 8.9e-6; %[Sec]
Tx_SRRC_Filter_length = 7; %SRRC Filter Length
N=length(r_A); % Length of received signal
Ts_A=1/Fs_A;   % # symbols, oversampling factor
time=Ts_A*N; t=Ts_A:Ts_A:time; % sampling interval & time vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1:                    Carrier Recovery                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Freq_Offset_A, theta_A] = Carrier_Recovery(r,Fs_A,Fc_A); %determine the input/output (Carrier is input, Theta/Freq_off is output)


% Spectrum Of Mystery Signal A

R_A = 20*log10(abs(fftshift(fft(r_A)))/length(r_A))+30 ; % spectrum of r_A
freq_ax = (-length(R_A)/2:length(R_A)/2-1)/(Ts_A*length(R_A) ) ; % frequency vector

figure(1)

subplot(2,1,1),plot(t,r_A,'r')
xlabel('Time [mSec]'),ylabel('Amplitude'),grid
axis([min(t) max(t) -10 10])
title('Mystery Signal A')

subplot(2,1,2),plot(freq_ax/1e3,R_A,'r')
xlabel('Frequency [KHz]'),ylabel('Amplitude [dBm]'),grid
axis([min(freq_ax)/1e3 max(freq_ax)/1e3 min(R_A) max(R_A)])
title('Power Spectrum Density of Mystery A Signal at IF frequency')

%%%%%

figure(2)

if (Fc_A - floor(Fc_A/Fs_A)*Fs_A) > (Fs_A/2)
    fz_A = Fs_A- (Fc_A - floor(Fc_A/Fs_A)*Fs_A); % freq. at receiver
else
    fz_A = (Fc_A - floor(Fc_A/Fs_A)*Fs_A); % freq. at receiver
end

subplot(2,1,1),plot(t,r_A,'r')
xlabel('Time [mSec]'),ylabel('Amplitude'),grid
axis([min(t) max(t) -10 10])
title('Mystery Signal A')

subplot(2,1,2),plot(freq_ax/1e3,R_A,'r')
xlabel('Frequency [KHz]'),ylabel('Amplitude [dBm]'),grid
axis([min(freq_ax)/1e3 max(freq_ax)/1e3 min(R_A) max(R_A)])
xlim([(fz_A-2e3)/1e3, (fz_A+2e3)/1e3])
title('Power Spectrum Density of Mystery A Signal at IF frequency (Zoom-In)')

%%%%%

figure(3)

plot(t,theta_A),
title('Phase Tracking via the Costas Loop')
xlabel('Time [Sec]'); ylabel('Phase Offset [Radian]'), grid

disp('The offset frequency of Mystery A Signal with respect to Fc is about [Hz]: ')
disp(Freq_Offset_A);
disp('>> 1- Carrier Recovery.');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2:                    DownConversion                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[MysteryA_Signal_BaseBand, MysteryA_Signal_BaseBand_raw_phase] = DownConversion(r_A, Fs_A, Fc_A, theta_A);
disp('>> 2- Downconversion.');

figure(4)

% Spectrum

R_BB_A = fft(MysteryA_Signal_BaseBand)/(length(MysteryA_Signal_BaseBand)) ; % spectrum of rsc
R_BB_A = 20*log10(abs(fftshift(R_BB_A)/time))+30;
freq_BB_ax = (-length(R_BB_A)/2:length(R_BB_A)/2-1)/(Ts_A*length(R_BB_A) ) ; % frequency vector

subplot(2,1,1),plot(t,MysteryA_Signal_BaseBand,'b')
xlabel('Time [mSec]'),ylabel('Amplitude'),grid
axis([min(t) max(t) -10 10])
title('Baseband Mystery Signal A')


subplot(2,1,2),plot(freq_BB_ax/1e3,R_BB_A,'b');
xlabel('Frequency [KHz]'),ylabel('Amplitude [dBm]'),grid;
axis([min(freq_BB_ax)/1e3 max(freq_BB_ax)/1e3 min(R_BB_A) max(R_BB_A)]);
title('Power Spectrum Density of Baseband Mystery A Signal');

%=================================> Constellation Diagram

%Calculate the Difference between the sampling frequency and the sample rate e.g multiply
%symbol period by the sample rate. This allows us to get the correct number
%of samples per symbols (Highest peak). hence the reason why we have to upsample and
%downsample. This was done in order not to have a wrong constellation

Fs_R_ratio = 369/50;
Upsampling_factor_A = 50;
DownSample_factor_A = 369;
Display_buffer_A = 50;


% This line  performs interpolation on the MysteryA baseband and the upsampling factor  
temp_BB_Ar = interp(MysteryA_Signal_BaseBand_raw_phase,Upsampling_factor_A);
DownSampled_IQ_Ar = temp_BB_Ar(1:DownSample_factor_A:end);


%==============> Constellation Diagram With Phase Traking
temp_BB_A = interp(MysteryA_Signal_BaseBand,Upsampling_factor_A);
DownSampled_IQ_A = temp_BB_A(1:DownSample_factor_A:end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3:                    Clock Recovery                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[Rx_Symbols_A, tau] = Clock_Recovery(MysteryA_Signal_BaseBand,beta_A,Symbol_Period_A,Fs_A,Tx_SRRC_Filter_length); %tau was added to check if the clock will work well
disp('>> 3- Clock Recovery.');%, pause(10)

figure(5)

subplot(2,1,1), plot(Rx_Symbols_A(1:end-10),'b.')    % plot constellation diagram
title('Constellation diagram for Mystery A Signal');
ylabel('Estimated symbol values'),grid
subplot(2,1,2), plot(tau(1:end-10))               % plot trajectory of tau
ylabel('Offset Estimates'), xlabel('iterations'),grid


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4:        Demodulation, Equalization and Text Reconstruction       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[reconstructed_message_A, Symbol_A_Eq, Correlator_Output, Eq_error_Vec, Eq_taps] = Demodulator_4PAM(Rx_Symbols_A);

figure(6)

plot(Symbol_A_Eq,'b.')    % plot constellation diagram
title('Constellation diagram for Mystery A Signal after Equalization');
ylabel('Estimated symbol values')
xlabel('iterations'),grid

figure(7)
plot(Eq_error_Vec ,'-*k')
ylabel('Equalization Error')
xlabel('Iteration'),grid

figure(8)

plot(Correlator_Output,'-*k')
ylabel('Synchronization Correlator Output')
xlabel('Symbol index'),grid

figure(9)

freqz(Eq_taps)
title('Equalizer Frequency Response');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Mystery B Receiver%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('mysteryB.mat')
r_B = r.';


%==========> Signal parameters
M = 4; 
refC = pammod(0:M-1,M); % 4-PAM Reference Constellation
Fs_B = 950e3; %[Hz] Sampling freq.
Fc_B = 1.2e6; %[Hz] Carrier freq.
beta_B = 0.26; % rolloff parameter for srrc
Symbol_Period_B = 7.5e-6; %[Sec]
Tx_SRRC_Filter_length_B = 4; %SRRC Filter Length
N_B=length(r_B); 
Ts_B=1/Fs_B;   % # symbols, oversampling factor
time_B=Ts_B*N_B; t_B=Ts_B:Ts_B:time_B; % sampling interval & time vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1:                    Carrier Recovery                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Freq_Offset_B, theta_B] = Carrier_Recovery(r_B,Fs_B,Fc_B);


% Spectrum Of Mystery Signal B

R_B = 20*log10(abs(fftshift(fft(r_B)))/length(r_B))+30 ; % spectrum of r_B
freq_ax_B = (-length(R_B)/2:length(R_B)/2-1)/(Ts_B*length(R_B) ) ; % frequency vector

figure(10)

subplot(2,1,1),plot(t_B,r_B,'r')
xlabel('Time [mSec]'),ylabel('Amplitude'),grid
axis([min(t_B) max(t_B) -10 10])
title('Mystery Signal B')

subplot(2,1,2),plot(freq_ax_B/1e3,R_B,'r')
xlabel('Frequency [KHz]'),ylabel('Amplitude [dBm]'),grid
axis([min(freq_ax_B)/1e3 max(freq_ax_B)/1e3 min(R_B) max(R_B)])
title('Power Spectrum Density of Mystery B Signal at IF frequency')
%%%%%


figure(11)

if (Fc_B - floor(Fc_B/Fs_B)*Fs_B) > (Fs_B/2)
    fz_B = Fs_B - (Fc_B - floor(Fc_B/Fs_B)*Fs_B); % freq. at receiver
else
    fz_B = (Fc_B - floor(Fc_B/Fs_B)*Fs_B); % freq. at receiver
end

subplot(2,1,1),plot(t_B,r_B,'r')
xlabel('Time [mSec]'),ylabel('Amplitude'),grid
axis([min(t_B) max(t_B) -10 10])
title('Mystery Signal B')

subplot(2,1,2),plot(freq_ax_B/1e3,R_B,'r')
xlabel('Frequency [KHz]'),ylabel('Amplitude [dBm]'),grid
axis([min(freq_ax_B)/1e3 max(freq_ax_B)/1e3 min(R_B) max(R_B)])
xlim([(fz_B-2e3)/1e3, (fz_B+2e3)/1e3])
title('Power Spectrum Density of Mystery B Signal at IF frequency')
%%%%%

figure(12)

plot(t_B,theta_B),
title('Phase Tracking via the Costas Loop')
xlabel('Time [Sec]'); ylabel('Phase Offset [Radian]'), grid

disp('The offset frequency with respect to true Fc is about [Hz]: ')
disp(Freq_Offset_B);
disp('>> 1- Carrier Recovery.');%, pause(10)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2:                    DownConversion                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[MysteryB_Signal_BaseBand, MysteryB_Signal_BaseBand_raw_phase] = DownConversion(r_B, Fs_B, Fc_B, theta_B);
disp('>> 2- Downconversion.');%, pause(10)

figure(13)

% Spectrum

R_BB_B = fft(MysteryB_Signal_BaseBand)/(length(MysteryB_Signal_BaseBand)) ; % spectrum of rsc
R_BB_B = 20*log10(abs(fftshift(R_BB_B)/time_B))+30;
freq_BB_ax_B = (-length(R_BB_B)/2:length(R_BB_B)/2-1)/(Ts_B*length(R_BB_B) ) ; % frequency vector

subplot(2,1,1),plot(t_B,MysteryB_Signal_BaseBand,'b')
xlabel('Time [mSec]'),ylabel('Amplitude'),grid
axis([min(t_B) max(t_B) -10 10])
title('Baseband Mystery Signal B')


subplot(2,1,2),plot(freq_BB_ax_B/1e3,R_BB_B,'b');
xlabel('Frequency [KHz]'),ylabel('Amplitude [dBm]'),grid;
axis([min(freq_BB_ax_B)/1e3 max(freq_BB_ax_B)/1e3 min(R_BB_B) max(R_BB_B)]);
title('Power Spectrum Density of Baseband Mystery B Signal');

%=================================> Constellation Diagram
Fs_R_ratio = 1241/250;
Upsampling_factor_B = 250;
DownSample_factor_B = 1241;
Display_buffer_B = 50;

%==============> Constellation Diagram Without Phase Traking

temp_BB_Br = interp(MysteryB_Signal_BaseBand_raw_phase,Upsampling_factor_B);
DownSampled_IQ_Br = temp_BB_Br(1:DownSample_factor_B:end);


%==============> Constellation Diagram With Phase Traking

temp_BB_B = interp(MysteryB_Signal_BaseBand,Upsampling_factor_B);
DownSampled_IQ_B = temp_BB_B(1:DownSample_factor_B:end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3:                    Clock Recovery                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[Rx_Symbols_B, tau] = Clock_Recovery(MysteryB_Signal_BaseBand,beta_B,Symbol_Period_B,Fs_B,Tx_SRRC_Filter_length_B);
disp('>> 3- Clock Recovery.');%, pause(10)

figure(14)

subplot(2,1,1), plot(Rx_Symbols_B(1:end-10),'b.')    % plot constellation diagram
title('Constellation diagram for Mystery B Signal');
ylabel('Estimated symbol values'),grid
subplot(2,1,2), plot(tau(1:end-10))               % plot trajectory of tau
ylabel('Offset Estimates'), xlabel('iterations'),grid


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4:        Demodulation, Equalization and Text Reconstruction       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[reconstructed_message_B, Symbol_B_Eq, Correlator_Output_B, Eq_error_Vec_B, Eq_taps] = Demodulator_4PAM(Rx_Symbols_B);

figure(15)

plot(Symbol_B_Eq,'b.')    % plot constellation diagram
title('Constellation diagram for Mystery B Signal after Equalization');
ylabel('Estimated symbol values')
xlabel('iterations'),grid


figure(16)
plot(Eq_error_Vec_B ,'-*r')
ylabel('Equalization Error')
xlabel('Iterations'),grid

figure(17)

plot(Correlator_Output_B,'-*r')
ylabel('Synchronization Correlator Output')
xlabel('Iterations'),grid

figure(18)

freqz(Eq_taps)
title('Equalizer Frequency Response');
%##########################################################################


%%%%%%%%%%%%%%%%%%%%%
% Mystery C Receiver%
%%%%%%%%%%%%%%%%%%%%%

% clear all;
% close all;
% clc;
load('mysteryC.mat')


r_C = r.';


%==========> Signal parameters
M = 4; 
refC = pammod(0:M-1,M); % 4-PAM Reference Constellation
Fs_C = 819e3; %[Hz] Sampling freq.
Fc_C = 2.2e6; %[Hz] Carrier freq.
beta_C = 0.32; % rolloff parameter for srrc
Symbol_Period_C = 8.14e-6; %[Sec]
Tx_SRRC_Filter_length = 6; %SRRC Filter Length
N=length(r_C); 
Ts_C=1/Fs_C;   % # symbols, oversampling factor
time=Ts_C*N; t=Ts_C:Ts_C:time; % sampling interval & time vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1:                    Carrier Recovery                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Freq_Offset_C, theta_C] = Carrier_Recovery(r,Fs_C,Fc_C);


% Spectrum Of Mystery Signal C

R_C = 20*log10(abs(fftshift(fft(r_C)))/length(r_C))+30 ; % spectrum of r_C
freq_ax = (-length(R_C)/2:length(R_C)/2-1)/(Ts_C*length(R_C) ) ; % frequency vector

figure(19)

subplot(2,1,1),plot(t,r_C,'r')
xlabel('Time [mSec]'),ylabel('Amplitude'),grid
axis([min(t) max(t) -10 10])
title('Mystery Signal C')

subplot(2,1,2),plot(freq_ax/1e3,R_C,'r')
xlabel('Frequency [KHz]'),ylabel('Amplitude [dBm]'),grid
axis([min(freq_ax)/1e3 max(freq_ax)/1e3 min(R_C) max(R_C)])
title('Power Spectrum Density of Mystery C Signal at IF frequency')
%%%%%

figure(20)

if (Fc_C - floor(Fc_C/Fs_C)*Fs_C) > (Fs_C/2)
    fz_C = Fs_C- (Fc_C - floor(Fc_C/Fs_C)*Fs_C); % freq. at receiver
else
    fz_C = (Fc_C - floor(Fc_C/Fs_C)*Fs_C); % freq. at receiver
end

subplot(2,1,1),plot(t,r_C,'r')
xlabel('Time [mSec]'),ylabel('Amplitude'),grid
axis([min(t) max(t) -10 10])
title('Mystery Signal C')

subplot(2,1,2),plot(freq_ax/1e3,R_C,'r')
xlabel('Frequency [KHz]'),ylabel('Amplitude [dBm]'),grid
axis([min(freq_ax)/1e3 max(freq_ax)/1e3 min(R_C) max(R_C)])
xlim([(fz_C-2e3)/1e3, (fz_C+2e3)/1e3])
title('Power Spectrum Density of Mystery C Signal at IF frequency')

%%%%%

figure(21)

plot(t,theta_C),
title('Phase Tracking via the Costas Loop')
xlabel('Time [Sec]'); ylabel('Phase Offset [Radian]'), grid

disp('The offset frequency of Mystery C Signal with respect to true Fc is about [Hz]: ')
disp(Freq_Offset_C);
disp('>> 1- Carrier Recovery.');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2:                    DownConversion                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[MysteryC_Signal_BaseBand, MysteryC_Signal_BaseBand_raw_phase] = DownConversion(r_C, Fs_C, Fc_C, theta_C);
disp('>> 2- Downconversion.');%, pause(10)

figure(22)

% Spectrum

R_BB_C = fft(MysteryC_Signal_BaseBand)/(length(MysteryC_Signal_BaseBand)) ; % spectrum of rsc
R_BB_C = 20*log10(abs(fftshift(R_BB_C)/time))+30;
freq_BB_ax = (-length(R_BB_C)/2:length(R_BB_C)/2-1)/(Ts_C*length(R_BB_C) ) ; % frequency vector

subplot(2,1,1),plot(t,MysteryC_Signal_BaseBand,'b')
xlabel('Time [mSec]'),ylabel('Amplitude'),grid
axis([min(t) max(t) -10 10])
title('Baseband Mystery Signal C')


subplot(2,1,2),plot(freq_BB_ax/1e3,R_BB_C,'b');
xlabel('Frequency [KHz]'),ylabel('Amplitude [dBm]'),grid;
axis([min(freq_BB_ax)/1e3 max(freq_BB_ax)/1e3 min(R_BB_C) max(R_BB_C)]);
title('Power Spectrum Density of Baseband Mystery C Signal');

%=================================> Constellation Diagram
Fs_R_ratio = 779/125;
Upsampling_factor_C = 125;
DownSample_factor_C = 779;
Display_buffer_C = 50;

%==============> Constellation Diagram Without Phase Traking
temp_BB_Cr = interp(MysteryC_Signal_BaseBand_raw_phase,Upsampling_factor_C);
DownSampled_IQ_Cr = temp_BB_Cr(1:DownSample_factor_C:end);

%==============> Constellation Diagram With Phase Traking
temp_BB_C = interp(MysteryC_Signal_BaseBand,Upsampling_factor_C);
DownSampled_IQ_C = temp_BB_C(1:DownSample_factor_C:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3:                    Clock Recovery                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[Rx_Symbols_C, tau] = Clock_Recovery(MysteryC_Signal_BaseBand,beta_C,Symbol_Period_C,Fs_C,Tx_SRRC_Filter_length);
disp('>> 3- Clock Recovery.');%, pause(10)


figure(23)

subplot(2,1,1), plot(Rx_Symbols_C(1:end-10),'b.')    % plot constellation diagram
title('Constellation diagram for Mystery C Signal');
ylabel('Estimated symbol values'),grid
subplot(2,1,2), plot(tau(1:end-10))               % plot trajectory of tau
ylabel('Offset Estimates'), xlabel('iterations'),grid


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4:        Demodulation, Equalization and Text Reconstruction       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[reconstructed_message_C, Symbol_C_Eq, Correlator_Output, Eq_error_Vec, Eq_taps] = Demodulator_4PAM(Rx_Symbols_C);

figure(24)

plot(Symbol_C_Eq,'b.')    % plot constellation diagram
title('Constellation diagram for Mystery C Signal after Equalization');
ylabel('Estimated symbol values')
xlabel('iterations'),grid


figure(25)
plot(Eq_error_Vec ,'-*k')
ylabel('Equalization Error')
xlabel('Symbols index'),grid

figure(26)

plot(Correlator_Output,'-*k')
ylabel('Synchronization Correlator Output')
xlabel('Symbol index'),grid

figure(27)

freqz(Eq_taps)
title('Equalizer Frequency Response');
%##########################################################################
disp('>> ================================================== <<');
disp('>> The reconstructed Message from Mystery A Signal is <<');
disp('>> ================================================== <<');
disp(reconstructed_message_A);
disp('>> ######################################################################## <<');


disp('>>                                                    <<');
disp('>> ================================================== <<');
disp('>> The reconstructed Message from Mystery B Signal is <<');
disp('>> ================================================== <<');
disp(reconstructed_message_B);
disp('>> ######################################################################## <<');


disp('>>                                                    <<');
disp('>> ================================================== <<');
disp('>> The reconstructed Message from Mystery C Signal is <<');
disp('>> ================================================== <<');

disp(reconstructed_message_C);
disp('>> ####################################################################### <<');
