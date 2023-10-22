bitstream = [1,0,1,1,1,1,0,0];

%BPSK. We assume M = 2, and this means that we have two symbols. 
% Since we have 8 symbols in our bitstream, duration of our signal must be
% 8Ts

frequency = 3; %Frequency of sinusoidal waves

ts = 1; % Symbol duration

sampling_rate = frequency * 20; %I took 10 times of Nyquist Rate as sampling rate.


t = linspace(0,8*ts,sampling_rate*8*ts);
carrier = exp((1i)*2*pi*frequency * t); %Our carrier. I defined it as a complex signal. I will use real part of it.

multiplier_bpsk = sqrt(2);
g = ones(1,ts*sampling_rate); 

modulated_signal = zeros(1,sampling_rate*8*ts);
constellation_bpsk = zeros(1,8);
for j=1:8
    if (bitstream(j) == 0)
        modulated_signal(1,((j-1)*ts*sampling_rate+1):j*ts*sampling_rate) = (-1) * multiplier_bpsk * g .* carrier(1,((j-1)*ts*sampling_rate+1):j*ts*sampling_rate); %pi degree shift is simply multiplting with -1 and I multiply with sqrt(2) to make average energy 1.
        constellation_bpsk(j) = sqrt(2/2)*(-1); %In order to create an orthonormal basis. 1 comes from sqrt(Eg/2). We know that Eg = 2.
    else 
        modulated_signal(1,((j-1)*ts*sampling_rate+1):j*ts*sampling_rate) = multiplier_bpsk * g .* carrier(1,((j-1)*ts*sampling_rate+1):j*ts*sampling_rate);
        constellation_bpsk(j) = sqrt(2/2);
    end 
end 

bpsk_modulated_pulse_stream = real(modulated_signal); %Our BPSK modulated signal.

%QPSK
% Now we have M = 4, it means that every 2 bit represents a symbol.

bitstream_qpsk = reshape(bitstream,[2,4]); %I am reshaping the array here
bitstream_qpsk = bitstream_qpsk';

t_qpsk = linspace(0,4*ts,sampling_rate*4*ts);

carrier_qpsk = exp((1i)*2*pi*frequency * t_qpsk);

qpsk_modulated_signal = zeros(1,sampling_rate*4*ts);

multiplier_qpsk = sqrt(2);

constellation_qpsk = zeros(1,4);
for k = 1:4
    if bitstream_qpsk(k,:) == [0,0]
         qpsk_modulated_signal(1,((k-1)*ts*sampling_rate+1):k*ts*sampling_rate) = multiplier_qpsk *  g .* carrier_qpsk(1,((k-1)*ts*sampling_rate+1):k*ts*sampling_rate); 
         constellation_qpsk(k) = sqrt(2/2);
    elseif bitstream_qpsk(k,:) == [0,1]
        qpsk_modulated_signal(1,((k-1)*ts*sampling_rate+1):k*ts*sampling_rate) = multiplier_qpsk *  g .* carrier_qpsk(1,((k-1)*ts*sampling_rate+1):k*ts*sampling_rate) * 1i; 
        constellation_qpsk(k) = sqrt(2/2)*1i;
    elseif bitstream_qpsk(k,:) == [1,0]
        qpsk_modulated_signal(1,((k-1)*ts*sampling_rate+1):k*ts*sampling_rate) = multiplier_qpsk * g .* carrier_qpsk(1,((k-1)*ts*sampling_rate+1):k*ts*sampling_rate) * (-1i); 
        constellation_qpsk(k) = sqrt(2/2)*(-1i);
    elseif bitstream_qpsk(k,:) == [1,1]
        qpsk_modulated_signal(1,((k-1)*ts*sampling_rate+1):k*ts*sampling_rate) = multiplier_qpsk * g .* carrier_qpsk(1,((k-1)*ts*sampling_rate+1):k*ts*sampling_rate) * (-1);
        constellation_qpsk(k) = sqrt(2/2)*(-1);
    end
end

qpsk_modulated_pulse_stream = real(qpsk_modulated_signal); %Our QPSK modulated signal

%4-PAM

% We have 4 symbols, which states that M = 4. And also we know that Eg = 1.
% To make average symbol energy unity, we need to take delta = sqrt(2/5) 

bitstream_4_pam = reshape(bitstream,[2,4]); %I am reshaping the array here
bitstream_4_pam = bitstream_4_pam';

t_4_pam = linspace(0,4*ts,sampling_rate*4*ts); 

carrier_4_pam = exp((1i)*2*pi*frequency * t_4_pam);

four_pam_signal = zeros(1,sampling_rate*4*ts); 

delta = sqrt(2/5);

constellation_4_pam = zeros(1,4);
for i= 1:4  
    if bitstream_4_pam(i,:) == [0,0]
         four_pam_signal(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate) = (-3)* delta * g .* carrier_4_pam(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate); 
         constellation_4_pam(i) = sqrt(1/2) * delta * (-3);
    elseif bitstream_qpsk(i,:) == [0,1]
        four_pam_signal(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate) = (-1) * delta* g .* carrier_4_pam(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate); 
        constellation_4_pam(i) = sqrt(1/2) * delta * (-1);
    elseif bitstream_qpsk(i,:) == [1,0]
        four_pam_signal(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate) = 3 * delta * g .* carrier_4_pam(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate); 
        constellation_4_pam(i) = sqrt(1/2) * delta * (3);
    elseif bitstream_qpsk(i,:) == [1,1]
        four_pam_signal(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate) =  1 * delta * g .* carrier_4_pam(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate);
        constellation_4_pam(i) = sqrt(1/2) * delta * (1);
    end
end

four_pam_signal = real(four_pam_signal); %Our 4-PAM Modulated Signal.

%16-QAM 
% With QAM I can send two symbols in one symbol through using
% orthogonality. So I need two symbols to transmit bitstream.

bitstream_16_qam= reshape(bitstream,[4,2]); %I am reshaping the array here
bitstream_16_qam = bitstream_16_qam';

t_16_qam = linspace(0,2*ts,sampling_rate*ts*2);

carrier_16_qam = exp((1i)*2*pi*frequency * t_16_qam); %Since I use two orthogonal sinusoidal waves, I can utilize imaginary part of the euler's identity here. 

delta_qam = sqrt(1/5);

qam_signal = zeros(1,sampling_rate*2*ts);

constellation_16_qam = zeros(1,2);
for i= 1:2  
    if bitstream_16_qam(i,:) == [0,0,0,0]
        qam_signal(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate) = (-3) * delta_qam * g .* real(carrier_16_qam(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate)) - (-3) * delta_qam * g .* imag(carrier_16_qam(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate)) ; 
        constellation_16_qam(i) = sqrt(1/2)*(-3)*delta_qam + sqrt(1/2)* delta_qam * (-3)* (1i);
    elseif bitstream_16_qam(i,:) == [0,0,0,1]
        qam_signal(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate) = (-3) * delta_qam * g .* real(carrier_16_qam(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate)) - (-1) * delta_qam * g .* imag(carrier_16_qam(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate)); 
        constellation_16_qam(i) = sqrt(1/2)*(-3)*delta_qam + sqrt(1/2)* delta_qam * (-1)* (1i);
    elseif bitstream_16_qam(i,:) == [0,0,1,0]
        qam_signal(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate) = (-3) * delta_qam * g .* real(carrier_16_qam(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate)) - 3 * delta_qam * g .* imag(carrier_16_qam(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate)); 
        constellation_16_qam(i) = sqrt(1/2)*(-3)*delta_qam + sqrt(1/2)* delta_qam * (3)* (1i);
    elseif bitstream_16_qam(i,:) == [0,0,1,1]
        qam_signal(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate) = (-3) * delta_qam * g .* real(carrier_16_qam(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate)) - delta_qam * g .* imag(carrier_16_qam(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate)); 
        constellation_16_qam(i) = sqrt(1/2)*(-3)*delta_qam + sqrt(1/2)* delta_qam * (1)* (1i);
    elseif bitstream_16_qam(i,:) == [0,1,0,0]
        qam_signal(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate) = (-1) * delta_qam * g .* real(carrier_16_qam(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate)) - (-3) * delta_qam * g .* imag(carrier_16_qam(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate)); 
        constellation_16_qam(i) = sqrt(1/2)*(-1)*delta_qam + sqrt(1/2)* delta_qam * (-3)* (1i);
    elseif bitstream_16_qam(i,:) == [0,1,0,1]
        qam_signal(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate) = (-1) * delta_qam * g .* real(carrier_16_qam(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate)) - (-1) * delta_qam * g .* imag(carrier_16_qam(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate)); 
        constellation_16_qam(i) = sqrt(1/2)*(-1)*delta_qam + sqrt(1/2)* delta_qam * (-1)* (1i);
    elseif bitstream_16_qam(i,:) == [0,1,1,0]
        qam_signal(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate) = (-1) * delta_qam * g .* real(carrier_16_qam(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate)) - 3 * delta_qam * g .* imag(carrier_16_qam(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate));
        constellation_16_qam(i) = sqrt(1/2)*(-1)*delta_qam + sqrt(1/2)* delta_qam * (3)* (1i);
    elseif bitstream_16_qam(i,:) == [0,1,1,1]
        qam_signal(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate) = (-1) * delta_qam * g .* real(carrier_16_qam(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate)) - 1 *  delta_qam * g .* imag(carrier_16_qam(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate)); 
        constellation_16_qam(i) = sqrt(1/2)*(-1)*delta_qam + sqrt(1/2)* delta_qam * (1)* (1i);
    elseif bitstream_16_qam(i,:) == [1,0,0,0]
        qam_signal(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate) = 3 * delta_qam * g .* real(carrier_16_qam(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate)) - (-3) * delta_qam * g .* imag(carrier_16_qam(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate)); 
        constellation_16_qam(i) = sqrt(1/2)*(3)*delta_qam + sqrt(1/2)* delta_qam * (-3)* (1i);
    elseif bitstream_16_qam(i,:) == [1,0,0,1]
        qam_signal(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate) = 3 * delta_qam * g .* real(carrier_16_qam(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate)) - (-1) * delta_qam * g .* imag(carrier_16_qam(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate));
        constellation_16_qam(i) = sqrt(1/2)*(3)*delta_qam + sqrt(1/2)* delta_qam * (-1)* (1i);
    elseif bitstream_16_qam(i,:) == [1,0,1,0]
        qam_signal(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate) = 3 * delta_qam *  g .* real(carrier_16_qam(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate)) - 3 * delta_qam * g .* imag(carrier_16_qam(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate)); 
        constellation_16_qam(i) = sqrt(1/2)*(3)*delta_qam + sqrt(1/2)* delta_qam * (3)* (1i);
    elseif bitstream_16_qam(i,:) == [1,0,1,1]
        qam_signal(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate) = 3 * delta_qam * g .* real(carrier_16_qam(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate)) - 1 *  delta_qam * g .* imag(carrier_16_qam(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate)); 
        constellation_16_qam(i) = sqrt(1/2)*(3)*delta_qam + sqrt(1/2)* delta_qam * (1)* (1i);
    elseif bitstream_16_qam(i,:) == [1,1,0,0]
        qam_signal(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate) =  delta_qam * g .* real(carrier_16_qam(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate)) - (-3) * delta_qam * g .* imag(carrier_16_qam(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate)); 
        constellation_16_qam(i) = sqrt(1/2)*(1)*delta_qam + sqrt(1/2)* delta_qam * (-3)* (1i);
    elseif bitstream_16_qam(i,:) == [1,1,0,1]
        qam_signal(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate) = delta_qam * g .* real(carrier_16_qam(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate)) - (-1) * delta_qam * g .* imag(carrier_16_qam(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate)); 
        constellation_16_qam(i) = sqrt(1/2)*(1)*delta_qam + sqrt(1/2)* delta_qam * (-1)* (1i);
    elseif bitstream_16_qam(i,:) == [1,1,1,0]
        qam_signal(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate) = delta_qam * g .* real(carrier_16_qam(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate)) - 3 * delta_qam * g .* imag(carrier_16_qam(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate)); 
        constellation_16_qam(i) = sqrt(1/2)*(1)*delta_qam + sqrt(1/2)* delta_qam * (3)* (1i);
    elseif bitstream_16_qam(i,:) == [1,1,1,1]
        qam_signal(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate) =  delta_qam * g .* real(carrier_16_qam(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate)) - 1 *  delta_qam *  g .* imag(carrier_16_qam(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate));
        constellation_16_qam(i) = sqrt(1/2)*(1)*delta_qam + sqrt(1/2)* delta_qam * (1)* (1i);
        
    end
end

%BFSK
%We have 2 symbols. To transmit bitstream we need 8 Ts. 
bitstream_bfsk = bitstream;
t_bfsk = linspace(0,8*ts,sampling_rate*8*ts);


bfsk_modulated_signal = zeros(1,sampling_rate*8*ts); 

multiplier_bfsk = sqrt(2); %To make average energy of constellation unity.

delta_f = 1/(2*ts);
constellation_bfsk = zeros(1,2);
for i = 1:8
    if (bitstream_bfsk(i) == 0)
        bfsk_modulated_signal(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate) = multiplier_bfsk * g .* exp((1i)*2*pi*frequency * t_bfsk(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate)) .* exp(1i*2*pi* delta_f .* t_bfsk(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate));
        constellation_bfsk(i) = 1i;
    else 
        bfsk_modulated_signal(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate) = multiplier_bfsk * g .* exp((1i)*2*pi*frequency * t_bfsk(1,((i-1)*ts*sampling_rate+1):i*ts*sampling_rate));
        constellation_bfsk(i) = 1;
    end 
end

bfsk_modulated_signal = real(bfsk_modulated_signal); %Our BFSK modulated signal.
subplot(5,1,1);
plot(t,bpsk_modulated_pulse_stream);
title('BPSK','interpreter','latex');
xlabel('Time','interpreter','latex');
ylabel("Amplitude",'interpreter','latex');

subplot(5,1,2);
plot(t_qpsk,qpsk_modulated_pulse_stream);
title('QPSK','interpreter','latex');
xlabel('Time','interpreter','latex');
ylabel("Amplitude",'interpreter','latex');

subplot(5,1,3);
plot(t_4_pam,four_pam_signal);
title('4-PAM','interpreter','latex');
xlabel('Time','interpreter','latex');
ylabel("Amplitude",'interpreter','latex');

subplot(5,1,4);
plot(t_16_qam, qam_signal);
title('16-QAM','interpreter','latex');
xlabel('Time','interpreter','latex');
ylabel("Amplitude",'interpreter','latex'); 

subplot(5,1,5);
plot(t_bfsk,bfsk_modulated_signal);
title('BFSK','interpreter','latex');
xlabel('Time','interpreter','latex');
ylabel("Amplitude",'interpreter','latex');

scatterplot(constellation_bpsk);
title('Constellation of BPSK','interpreter','latex');
xlabel('$\sqrt{\frac{2}{E_g}}g(t)cos(2\pi f_c t)$','interpreter','latex');
ylabel("$\sqrt{\frac{2}{E_g}}g(t)sin(2\pi f_c t)$",'interpreter','latex');


scatterplot(constellation_qpsk);
title('Constellation of QPSK','interpreter','latex');
xlabel('$\sqrt{\frac{2}{E_g}}g(t)cos(2\pi f_c t)$','interpreter','latex');
ylabel("$\sqrt{\frac{2}{E_g}}g(t)sin(2\pi f_c t)$",'interpreter','latex');

scatterplot(constellation_4_pam);
title('Constellation of 4-PAM','interpreter','latex');
xlabel('$\sqrt{\frac{E_g}{2}}\Delta$','interpreter','latex');


scatterplot(constellation_16_qam);
title('Constellation of 16-QAM','interpreter','latex');
xlabel('$\sqrt{\frac{2}{E_g}}g(t)cos(2\pi f_c t)$','interpreter','latex');
ylabel("$\sqrt{\frac{2}{E_g}}g(t)sin(2\pi f_c t)$",'interpreter','latex');

scatterplot(constellation_bfsk);
title('Constellation of BFSK','interpreter','latex');
xlabel('$\sqrt{\frac{E_g}{2}}\phi_0$','interpreter','latex');
ylabel("$\sqrt{\frac{E_g}{2}}\phi_1$",'interpreter','latex');
