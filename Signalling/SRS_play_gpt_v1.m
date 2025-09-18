% Clear previous figures and variables
clear;
close all;
clc
%%
% Define carrier and SRS configurations
carrier = nrCarrierConfig;
carrier.NSizeGrid = 273; % Number of Resource Blocks; 273 = 100MHz
carrier.SubcarrierSpacing = 30; % Subcarrier spacing in kHz
carrier.NSlot = 0; % Slot number for this example. 

carrier
%%
% Configure the SRS
srs = nrSRSConfig;
srs.NumSRSPorts = 1;            % Number of SRS antenna ports (1,2,4). 1 for UL Estimation for 1T4R typical UE
srs.SymbolStart = 6;            % Starting OFDM symbol within a slot. 'The last 6 symbols of the slot, or at any symbol location within the slot if resourceMapping-r16 is provided subject to UE capability'
srs.NumSRSSymbols = 1;          % Number of OFDM symbols allocated per slot (1,2,4). Up to 8 for SRS-PosResourceSet-r16
srs.ResourceType = 'periodic';  % Resource type ('periodic', 'semi-persistent','aperiodic'). Use 'aperiodic' to disable inter-slot frequency hopping
srs.SRSPeriod = [80 0];         % Periodicity and offset in slots. SRSPeriod(2) must be < SRSPeriod(1)
srs.FrequencyStart = 0;         % Frequency position of the SRS in BWP in RBs
srs.NRRC = 0;                   % Additional offset from FreqStart specified in blocks of 4 PRBs (0...67)
srs.CSRS = 63%5;%11;            % Bandwidth configuration C_SRS (0...63). It controls the allocated bandwidth to the SRS. Table 6.4.1.4.3-1: SRS bandwidth configuration.
srs.BSRS = 0;                   % Bandwidth configuration B_SRS (0...3). It controls the allocated bandwidth to the SRS
srs.BHop = 0;                   % Frequency hopping configuration (0...3). Set BHop < BSRS to enable frequency hopping
srs.KTC = 2%8%4;                % Comb number (2,4). Frequency density in subcarriers. 2 - every second subcarrier.

%srs.SRSPositioning = 1;

srs
%%
% Generate SRS indices and symbols
[srsIndices, srsInfo] = nrSRSIndices(carrier, srs);
srsSymbols = nrSRS(carrier, srs);
size(srsIndices)
size(srsSymbols)

% Print sizes of generated indices and symbols
disp('Size of srsIndices [(12 subcarrier / 2 KTC) x 272 PRB = 1632]:');
disp(size(srsIndices));
disp('Size and type of srsSymbols [compex double]:');
disp(size(srsSymbols));
disp(class(srsSymbols));

%%
% Create and populate resource grid
K = carrier.NSizeGrid * 12; % Total number of subcarriers
L = 14; % Number of OFDM symbols per slot
P = srs.NumSRSPorts; % Number of ports
resourceGrid = complex(zeros(K, L, P));
resourceGrid(srsIndices) = srsSymbols;


%%
% OFDM Modulation Parameters
ofdmInfo = nrOFDMInfo(carrier, 'CarrierFrequency', 3.5e9);
txWaveform = nrOFDMModulate(carrier, resourceGrid);

%%
% Define the sampling rate based on OFDM parameters
Fs = ofdmInfo.SampleRate;

% Generate a time vector based on the length of the txWaveform
t = (0:length(txWaveform)-1) / Fs;

% Plot the real part of the txWaveform
figure;
plot(t, real(txWaveform));
xlabel('Time (s)');
ylabel('Amplitude');
title('Real Part of Transmitted Waveform in Time Domain');
grid on;

% Optional: plot the imaginary part or magnitude
figure;
plot(t, imag(txWaveform));
xlabel('Time (s)');
ylabel('Amplitude');
title('Imaginary Part of Transmitted Waveform in Time Domain');
grid on;

figure;
plot(t, abs(txWaveform));
xlabel('Time (s)');
ylabel('Amplitude');
title('Magnitude of Transmitted Waveform in Time Domain');
grid on;




%%
% Frequency domain analysis
N = length(txWaveform);
freqVector = (-N/2:N/2-1)*(Fs/N); % Frequency vector for FFT

% Perform FFT and shift zero frequency component to the center
txWaveformFFT = fftshift(fft(txWaveform));

% Plot the amplitude of the subcarriers
figure;
plot(freqVector, abs(txWaveformFFT));
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Amplitude of Subcarriers in Frequency Domain');
grid on;

% Plot the phase of the subcarriers
figure;
plot(freqVector, angle(txWaveformFFT));
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
title('Phase of Subcarriers in Frequency Domain');
grid on;


% Plot the spectrogram of the waveform
samplerate = ofdmInfo.SampleRate; % Use SampleRate from ofdmInfo
nfft = ofdmInfo.Nfft; % Use Nfft from ofdmInfo
figure;
spectrogram(txWaveform(:,1), ones(nfft, 1), 0, nfft, 'centered', samplerate, 'yaxis', 'MinThreshold', -130);
title('Spectrogram of 5G Uplink Baseband Waveform');

% 4. Plot the resource grid with SRS location
figure;
cmap = parula(64);

% Full RB - port 0
subplot(2,2,[1 3]);
hold on;
image(40*abs(resourceGrid(:,:,1)));
axis(gca,'xy');
colormap(cmap);
for i = 2:14
  line([i-0.5 i-0.5],[0 carrier.NSizeGrid*12],'Color','white');
end
xlim([0.5 14.5]);
ylim([0 carrier.NSizeGrid*12]);
set(gca,'xtick',0:14);
set(gca,'xticklabel',{'','0','1','2','3','4','5','6','7','8','9','10','11','12','13'});
set(gca,'ytick',[1 carrier.NSizeGrid*12]);
set(gca,'yticklabel',{'0',num2str(carrier.NSizeGrid-1)});
ylabel('RB');
tmpStr = sprintf('Full RB - port 0\nKTC=%d, CSRS=%d, BSRS=%d, BHop=%d', srs.KTC, srs.CSRS, srs.BSRS, srs.BHop);
title(tmpStr);
hold off;

% First RB - port 0
subplot(2,2,2);
hold on;
% Plot only the first RB for the first port
image(40*abs(resourceGrid(1:12,:,1))); 
axis(gca,'xy');
colormap(cmap);
for i = 2:14
  line([i-0.5 i-0.5],[0 12],'Color','white');
end
for j = 1:12
  line([0 15],[j+0.5 j+0.5],'Color','white');
end
xlim([0.5 14.5]);
ylim([0.5 12.5]);
set(gca,'xtick',0:14);
set(gca,'xticklabel',{'','0','1','2','3','4','5','6','7','8','9','10','11','12','13'});
set(gca,'ytick',0:12);
set(gca,'yticklabel',{'','0','1','2','3','4','5','6','7','8','9','10','11'});
ylabel('Subcarrier');
title('First RB - port 0');
hold off;


% Spectrogram of the waveform
figure;
spectrogram(txWaveform(:,1), ones(nfft, 1), 0, nfft, 'centered', samplerate, 'yaxis', 'MinThreshold', -130);
title('Spectrogram of 5G Uplink Baseband Waveform');

