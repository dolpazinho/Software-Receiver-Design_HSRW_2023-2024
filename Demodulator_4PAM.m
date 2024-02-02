%function [received_message] = Demodulator_4PAM(xs)
function [reconstructed_message, Symbol_A_Eq, Correlator_Output, Eq_error_Vec, Eq_taps] = Demodulator_4PAM(xs)

alphabet = [-3, -1, 1, 3]; %Initialise our symbol alphabets

%Parameters of the equalizer
Equalizer_order = 3; % Order must be odd number
Eq_delay = floor(Equalizer_order/2)+1;
Eq_taps = zeros(1,Equalizer_order);
Eq_error_Vec = [];
Symbol_A_hat = [];
Symbol_A_Eq = [];
mu = 0.01;

%===================>  Finding start of Message using preamble 
% (Detect the start of preamble)

Preamble_symbols = letters2pam('0x0 This is is the Frame Header 1y1');% Convert 
% the message into PAM symbol (letters2pam)

%======================================= Correlation ============================
xsc                    = filter(fliplr(Preamble_symbols),1,xs);
[max_value,istart]     = max(abs(xsc)); 
Correlator_Output = abs(xsc);

Length_Data_block = 500;  % This is to locate the beginning of our message. 
Length_Preamble = length(Preamble_symbols);
Length_full_frame = Length_Preamble + Length_Data_block;

% Identify the beginning of the message
First_preamble = istart-Length_Preamble + 1;
First_preamble = mod(First_preamble,Length_full_frame) + Length_full_frame - Eq_delay +1;

length_usable_xs = length(xs(First_preamble:end))-mod(length(xs(First_preamble:end)),Length_full_frame) ; % Starting point length

Number_frame = length_usable_xs/Length_full_frame;

xs = [xs, 0];
xs = xs(First_preamble:First_preamble+length_usable_xs-1+1);

%=============== Equalizer ===============================================%

for i_frame = 0:Number_frame-1 % we implemented our equalizer frame by frame

    % LMS For the preamble
    for i_Preamble = Eq_delay:Length_Preamble

        Equalizer_TDL = xs(((i_frame*Length_full_frame) + i_Preamble - Eq_delay + 1) : ((i_frame*Length_full_frame) + i_Preamble + Eq_delay - 1) );
        Eq_error = Preamble_symbols(i_Preamble) - (Equalizer_TDL * Eq_taps');
        Eq_taps = Eq_taps + mu * Eq_error * Equalizer_TDL;
        Eq_error_Vec = [Eq_error_Vec, Eq_error];

    end
    
    % Decision Directed for the Data
    for i_Symbol = Eq_delay-1:Length_Data_block-Eq_delay + 2

        Equalizer_TDL = [xs((i_frame*Length_full_frame) + Length_Preamble + i_Symbol  - Eq_delay + 1 : (i_frame*Length_full_frame) + Length_Preamble + i_Symbol + Eq_delay - 1)];
        Eq_Output = Equalizer_TDL * Eq_taps';
        Eq_error = quantalph(Eq_Output,alphabet) - Eq_Output;
        Eq_taps = Eq_taps + mu * Eq_error * Equalizer_TDL;
        Eq_error_Vec = [Eq_error_Vec, Eq_error];

        Symbol_A_Eq = [Symbol_A_Eq, Eq_Output];
        Symbol_A_hat = [Symbol_A_hat, quantalph(Eq_Output,alphabet)];
        
    end
end

%Reconstruct the message from the decision output
reconstructed_message = pam2letters(Symbol_A_hat);





