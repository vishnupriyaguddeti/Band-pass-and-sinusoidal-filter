close all
clear all 
fclose('all');

if (exist('__octave_config_info__')~=0) && (length(strfind(computer(),'linux'))>0)  % octave specific and linux specific
    graphics_toolkit('fltk')
end

% The size of the figures displayed and the font size in the figures were adjusted for a 4k display (3840 x 2160).
if (exist('__octave_config_info__')~=0) && (length(strfind(computer(),'w64'))>0)  || (length(strfind(computer(),'w32'))>0) % octave specific and windows specific
    figures_position_size=[200,300,1300,1000]; %[startx,starty,width,height]
    figures_font_size=30;
end

% load libraries
if exist('__octave_config_info__')~=0 % automatically determine if we are running Octave version 5.x or Matlab
    pkg load signal % Octave specific
end

if exist('__octave_config_info__')~=0 % octave specific
    randn('seed', 1234);
else
    rng('default')
    rng(1234);
end


%Reading input signal and calculating FFT, power spectrum
[Ip_Signal,Fs] = audioread('assign1_speech_with_bandpass_noise_and_tonal_noise.wav');
Ip_len = length(Ip_Signal);
noverlap=Fs/2;

FFT_IpSignal = fft(Ip_Signal);
FFT_IpSignal = abs(FFT_IpSignal);
FFT_IpSignal= FFT_IpSignal.*FFT_IpSignal/length(FFT_IpSignal);
if exist('__octave_config_info__')~=0 % octave specific
    [Pxx,W] = cpsd(Ip_Signal,Ip_Signal,hann(Fs),0,Fs,'twosided'); % cpsd() normalizes window rms value
else
    [Pxx,W] = cpsd(Ip_Signal,Ip_Signal,hann(Fs),noverlap,Fs,'twosided'); % cpsd() normalizes window rms value
end
Pxx=Pxx*2*pi;

freq=(0:1:Fs-1)/Fs;
n=0:1:Ip_len-1;

%Input Signal Plot
if (exist('__octave_config_info__')~=0) && (length(strfind(computer(),'w64'))>0)  || (length(strfind(computer(),'w32'))>0) % octave specific and windows specific
    figure(1,'position',figures_position_size); %[startx,starty,width,height]
else
    figure(1)
end
subplot(3,1,1)
plot(n, Ip_Signal)
if (exist('__octave_config_info__')~=0) && (length(strfind(computer(),'w64'))>0)  || (length(strfind(computer(),'w32'))>0) % octave specific and windows specific
    set(gca, 'fontsize', figures_font_size)
end
title('Input signal')
ylabel('x[n]')
xlabel('n')
grid

subplot(3,1,2)
plot(n, FFT_IpSignal)
if (exist('__octave_config_info__')~=0) && (length(strfind(computer(),'w64'))>0)  || (length(strfind(computer(),'w32'))>0) % octave specific and windows specific
    set(gca, 'fontsize', figures_font_size)
end
title('FFT |X[k]|^2/N')
ylabel('FFT(Input)')
xlabel('k')
grid
xlim auto
%ylim auto
ylim([0 1])

subplot(3,1,3)
plot(freq, 10*log10(Pxx))
if (exist('__octave_config_info__')~=0) && (length(strfind(computer(),'w64'))>0)  || (length(strfind(computer(),'w32'))>0) % octave specific and windows specific
    set(gca, 'fontsize', figures_font_size)
end
title('PSD Pxx[f]')
ylabel('Pxx')
xlabel('freq')
grid
xlim auto
ylim auto
%ylim([0 0.5])

%Stop band elliptical filter
[b, a] = ellip(8,1,100,[0.30 0.49],'stop');

[H, w] = freqz(b,a,512,'twosided');
[h, nh] = impz(b,a,Fs); % same length as PSD estim, for comparison with hmodel


if (exist('__octave_config_info__')~=0) && (length(strfind(computer(),'w64'))>0)  || (length(strfind(computer(),'w32'))>0) % octave specific and windows specific
    figure(2,'position',figures_position_size); %[startx,starty,width,height]
else
    figure(2)
end
subplot(2,1,1)
%plot(w, 20*np.log10(np.abs(H)))
plot(w, abs(H))
if (exist('__octave_config_info__')~=0) && (length(strfind(computer(),'w64'))>0)  || (length(strfind(computer(),'w32'))>0) % octave specific and windows specific
    set(gca, 'fontsize', figures_font_size)
end
title('LTI freq. response |H(w)| (Elliptical Filter)')
ylabel('H')
xlabel('w')
grid

subplot(2,1,2)
plot(nh, h)
if (exist('__octave_config_info__')~=0) && (length(strfind(computer(),'w64'))>0)  || (length(strfind(computer(),'w32'))>0) % octave specific and windows specific
    set(gca, 'fontsize', figures_font_size)
end
%stem(nh, h)
title('LTI impulse response h[n] (Elliptical Filter)')
ylabel('h')
xlabel('n')
grid

if (exist('__octave_config_info__')~=0) && (length(strfind(computer(),'w64'))>0)  || (length(strfind(computer(),'w32'))>0) % octave specific and windows specific
    figure(3,'position',figures_position_size); %[startx,starty,width,height]
else
    figure(3)
end
zplane(b,a)
if (exist('__octave_config_info__')~=0) && (length(strfind(computer(),'w64'))>0)  || (length(strfind(computer(),'w32'))>0) % octave specific and windows specific
    set(gca, 'fontsize', figures_font_size)
end
title('z-plane, poles and zeros of LTI H(z) (Elliptical Filter)')


%Input signal filtered from bandpass noise
[Noise_Filt, zf] = filter(b, a, Ip_Signal);

FFT_Ellip = fft(Noise_Filt);
FFT_Ellip = abs(FFT_Ellip);
FFT_Ellip= FFT_Ellip.*FFT_Ellip/length(FFT_Ellip);

if exist('__octave_config_info__')~=0 % octave specific
    [P_ellip,W] = cpsd(Noise_Filt,Noise_Filt,hann(Fs),0,Fs,'twosided'); % cpsd() normalizes window rms value
else
    [P_ellip,W] = cpsd(Noise_Filt,Noise_Filt,hann(Fs),noverlap,Fs,'twosided'); % cpsd() normalizes window rms value
end
P_ellip=P_ellip*2*pi; % Matlab cpsd values are 2pi lower than they should be

if (exist('__octave_config_info__')~=0) && (length(strfind(computer(),'w64'))>0)  || (length(strfind(computer(),'w32'))>0) % octave specific and windows specific
    figure(4,'position',figures_position_size); %[startx,starty,width,height]
else
    figure(4)
end
subplot(3,1,1)
plot(n, Noise_Filt)
if (exist('__octave_config_info__')~=0) && (length(strfind(computer(),'w64'))>0)  || (length(strfind(computer(),'w32'))>0) % octave specific and windows specific
    set(gca, 'fontsize', figures_font_size)
end
title('Ellip filtered x[n]')
ylabel('x1[n]')
xlabel('n')
grid
ylim([-0.8 0.8])

subplot(3,1,2)
plot(n, FFT_Ellip)
if (exist('__octave_config_info__')~=0) && (length(strfind(computer(),'w64'))>0)  || (length(strfind(computer(),'w32'))>0) % octave specific and windows specific
    set(gca, 'fontsize', figures_font_size)
end
title('Ellip filtered FFT |X1[k]|^2/N')
ylabel('FFT-X1')
xlabel('k')
grid
xlim auto
ylim([0 1*(10^-3)])
%ylim([0 (6*(10^-4))])

subplot(3,1,3)
plot(freq, 10*log10(P_ellip))
if (exist('__octave_config_info__')~=0) && (length(strfind(computer(),'w64'))>0)  || (length(strfind(computer(),'w32'))>0) % octave specific and windows specific
    set(gca, 'fontsize', figures_font_size)
end
title('PSD P ellip[f]')
ylabel('P ellip')
xlabel('freq')
grid
xlim auto
ylim auto


%Notch filters to remove sinusoidal noise
w0 = [0.227;(0.5034);(0.6268);(0.7652)];
r = [0.94;0.92;0.96;0.93];
%bw = [(0.015*2);(0.2);(0.2);(0.2)];
Prev_NotchFilt = Noise_Filt;
for i = 1:length(w0)
    c = [1 -2*cos(pi*w0(i)) 1];
    norm_num = 2-2*cos(pi*w0(i));
    c = c/norm_num ;
    d = [1 -2*r(i)*cos(pi*w0(i)) r(i)*r(i)];
    norm_den = 1-(2*r(i)*cos(pi*w0(i)))+(r(i)*r(i));
    d = d/norm_den;
    
    [H, w] = freqz(c,d,512,'twosided');
    [h, nh] = impz(c,d,Fs);
    
    if (exist('__octave_config_info__')~=0) && (length(strfind(computer(),'w64'))>0)  || (length(strfind(computer(),'w32'))>0) % octave specific and windows specific
        figure(4+i,'position',figures_position_size); %[startx,starty,width,height]
    else
        figure(4+i)
    end
    subplot(3,1,1)
    %plot(w, 20*np.log10(np.abs(H)))
    plot(w, abs(H))
    if (exist('__octave_config_info__')~=0) && (length(strfind(computer(),'w64'))>0)  || (length(strfind(computer(),'w32'))>0) % octave specific and windows specific
        set(gca, 'fontsize', figures_font_size)
    end
    title('LTI freq. response |H(w)| (Notch Filter)')
    ylabel('H')
    xlabel('w')
    grid
    
    subplot(3,1,2)
    plot(nh, h)
    if (exist('__octave_config_info__')~=0) && (length(strfind(computer(),'w64'))>0)  || (length(strfind(computer(),'w32'))>0) % octave specific and windows specific
        set(gca, 'fontsize', figures_font_size)
    end
    %stem(nh, h)
    title('LTI impulse response h[n] (Notch Filter)')
    ylabel('h')
    xlabel('n')
    grid
    
    subplot(3,1,3)
    zplane(c,d)
    if (exist('__octave_config_info__')~=0) && (length(strfind(computer(),'w64'))>0)  || (length(strfind(computer(),'w32'))>0) % octave specific and windows specific
        set(gca, 'fontsize', figures_font_size)
    end
    title('z-plane, poles and zeros of LTI H(z)')
    
    Notch_Filt = filter(c,d,Prev_NotchFilt);
    Prev_NotchFilt = Notch_Filt;
end

Op_Signal = Notch_Filt;

if (exist('__octave_config_info__')~=0) && (length(strfind(computer(),'w64'))>0)  || (length(strfind(computer(),'w32'))>0) % octave specific and windows specific
    figure(9,'position',figures_position_size); %[startx,starty,width,height]
else
    figure(9)
end
plot(n,Op_Signal);
if (exist('__octave_config_info__')~=0) && (length(strfind(computer(),'w64'))>0)  || (length(strfind(computer(),'w32'))>0) % octave specific and windows specific
    set(gca, 'fontsize', figures_font_size)
end
title('Output signal y1[n] with Transient')
ylabel('y1[n]')
xlabel('n')
grid


%Transient removal
Trans_Len = 950;
Op_Signal = Op_Signal(Trans_Len:end);
n = n(Trans_Len:end);

Audio_Signal = Op_Signal;
scaled = (Audio_Signal/max(abs(Audio_Signal)))*0.999;
audiowrite('Filtered_Audio_Signal_Assignment1.wav', scaled, Fs)

%Op_Signal = filter(b,a,Ip_Signal);
FFT_Op = fft(Audio_Signal);
FFT_Op = abs(FFT_Op);
FFT_Op= FFT_Op.*FFT_Op/length(FFT_Op);

if exist('__octave_config_info__')~=0 % octave specific
    [Pyy,W] = cpsd(Op_Signal,Op_Signal,hann(Fs),0,Fs,'twosided'); % cpsd() normalizes window rms value
else
    [Pyy,W] = cpsd(Op_Signal,Op_Signal,hann(Fs),noverlap,Fs,'twosided'); % cpsd() normalizes window rms value
end
Pyy=Pyy*2*pi; % Matlab cpsd values are 2pi lower than they should be

if (exist('__octave_config_info__')~=0) && (length(strfind(computer(),'w64'))>0)  || (length(strfind(computer(),'w32'))>0) % octave specific and windows specific
    figure(10,'position',figures_position_size); %[startx,starty,width,height]
else
    figure(10)
end
subplot(3,1,1)
plot(n, Audio_Signal)
if (exist('__octave_config_info__')~=0) && (length(strfind(computer(),'w64'))>0)  || (length(strfind(computer(),'w32'))>0) % octave specific and windows specific
    set(gca, 'fontsize', figures_font_size)
end
title('Output signal y[n]')
ylabel('y[n]')
xlabel('n')
grid

subplot(3,1,2)
plot(n, FFT_Op)
if (exist('__octave_config_info__')~=0) && (length(strfind(computer(),'w64'))>0)  || (length(strfind(computer(),'w32'))>0) % octave specific and windows specific
    set(gca, 'fontsize', figures_font_size)
end
title('FFT |Y[k]|^2/N')
ylabel('FFT[Output]-Y')
xlabel('k')
grid
ylim([0 0.0012])

subplot(3,1,3)
plot(freq, 10*log10(Pyy))
if (exist('__octave_config_info__')~=0) && (length(strfind(computer(),'w64'))>0)  || (length(strfind(computer(),'w32'))>0) % octave specific and windows specific
    set(gca, 'fontsize', figures_font_size)
end
title('PSD Pyy[f]')
ylabel('Pyy')
xlabel('freq')
grid
ylim([-160 -20])