function [Mystery_Signal_BaseBand, Mystery_Signal_BaseBand_raw_phase] = DownConversion(r, Fs, Fc, theta)


time = length(r)/Fs ; 
Ts=1/Fs; % sampling interval and time base
t= Ts : Ts : time ; %lent=length ( t ) ; % define a "time" vector

gamma=0; % interpolates the factorial function.
phi =theta; % freq & phase offset compensation

if (Fc - floor(Fc/Fs)*Fs) > (Fs/2)
    f0 = Fs- (Fc - floor(Fc/Fs)*Fs); % freq. at receiver
else
    f0 = (Fc - floor(Fc/Fs)*Fs); % freq. at receiver
end

% To check carrier with phase composition
c =cos (2* pi *( f0+gamma)* t+phi ) ; % create cosine for demod/ downconversion
x1= r .* c; % demod received signal
fbe=[0 0.2 0.4 1 ] ; fa=[1 1 0 0 ] ; fl = 100; % low pass filter design
b=remez ( fl , fbe , fa ) ; % impulse response of LPF
m1= filter (b ,1 , x1 ) ; % LPF the demodulated signal


%===============================> AGC used to adjust the amplititude output to a constant value, 
% even if the input voltage decreases or increases.
N=length(m1);
ds=0.6;                            % desired power of output / signal energy
mu=0.0001;                          % Algorithm step-size

a1=zeros(1,N); a1(1)=1;              % initialize AGC parameter
s1=zeros(1,N);
 
for k=1:N-1
  s1(k)=a1(k)*m1(k);                  % normalize by a to get s
  a1(k+1)=a1(k)-mu*(s1(k)^2-ds);      % adaptive update of a(k)
end

m1 = s1;


Mystery_Signal_BaseBand = m1;

%==========================================================================

% To check carrier without phase composition
cr =cos (2* pi *( f0+gamma)* t ) ; % create cosine for demod/ downconversion
x1r= r .* cr; % demod received s ignal
fbe=[0 0.1 0.3 1 ] ; fa=[1 1 0 0 ] ; fl = 100; % low pass filter design
b=remez ( fl , fbe , fa ) ; % impulse response of LPF
m1r= filter (b ,1 , x1r ) ; % LPF the demodulated signal


Mystery_Signal_BaseBand_raw_phase = m1r;

