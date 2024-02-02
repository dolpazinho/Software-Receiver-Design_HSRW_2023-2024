% Clock Recovery Using Maximizing Output Power

function [xs, tausave] = Clock_Recovery(Mystery_Signal_BaseBand,beta,Symbol_Period,Fs, l)


% calculating number of symbol in our mystery signal
m = Symbol_Period*Fs;
n = floor(length(Mystery_Signal_BaseBand)/m);

matchfilt=srrc(l,beta,m,0);        % filter = pulse shape
x=conv(real(Mystery_Signal_BaseBand(1:end)),matchfilt);               % convolve signal with matched filter

x = x(1:n*m); % Ensure filter consistency so that we can have the same number of filter with the transmitter

% run clock recovery algorithm
tnow=l*m+1; tau=0; xs=zeros(1,n);   % initialize variables
tausave=zeros(1,n); tausave(1)=tau; i=0;
mu=0.01;                            % algorithm stepsize
delta=1;                          % time for derivative
while tnow<length(x)-l*m            % run iteration
  i=i+1;
  clc
  xs(i)=interpsinc(x,tnow+tau,l,beta);   % interp at tnow+tau
  x_deltap=interpsinc(x,tnow+tau+delta,l,beta);  % value to right
  x_deltam=interpsinc(x,tnow+tau-delta,l,beta);  % value to left
  dx=x_deltap-x_deltam;             % numerical derivative
  
  tau=tau+mu*dx*xs(i);              % algorithm update
  tnow=tnow+m; tausave(i)=tau;      % save for plotting
  
end


