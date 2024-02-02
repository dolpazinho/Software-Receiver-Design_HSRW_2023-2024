
%The implementation of this script as described in the SDR book Figure
%10.9

function [Freq_Offset, theta] = Carrier_Recovery(r,Fs,Fc)

%=======================> Signal parameters
N=length(r.'); 
Ts=1/Fs;   % No of symbols, oversampling factor
time=Ts*N; t=Ts:Ts:time; % sampling interval & time vector

%=====> costasloop.m simulate costas loop with input from pulrecsig.m
%r=rsc;                                % rsc from pulrecsig.m
fl=100; ff=[0 .1 .4 1]; fa=[1 1 0 0];   % Filter design
h=firpm(fl,ff,fa);                   
mu=.01;                              % algorithm stepsize


if (Fc - floor(Fc/Fs)*Fs) > (Fs/2)  %calculate estimated value 
    f0 = Fs- (Fc - floor(Fc/Fs)*Fs); % freq. at receiver
else
    f0 = (Fc - floor(Fc/Fs)*Fs); % freq. at receiver
end

theta=zeros(1,length(t)); theta(1)=0; % estimate vector
zs=zeros(1,fl+1); zc=zeros(1,fl+1);   % buffers for LPFs

for k=1:length(t)-1                  
    zs=[zs(2:fl+1), 2*r(k)*sin(2*pi*f0*t(k)+theta(k))];
    zc=[zc(2:fl+1), 2*r(k)*cos(2*pi*f0*t(k)+theta(k))];
    lpfs=fliplr(h)*zs'; lpfc=fliplr(h)*zc'; % output of filters
    theta(k+1)=theta(k)-mu*lpfs*lpfc;   % algorithm update
end

Freq_Offset = mean(diff(theta))/Ts/2/pi; % Estimate the freq Offset


