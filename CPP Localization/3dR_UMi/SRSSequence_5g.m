function [txFreqCoefficients, txWaveForm] = SRSSequence_5g(numSrsSymbols, symbolStartIdx, cSrs, combNumber)

%%--------------------------------------------------------------------------
% Author: Darius Chmielikaus (with contributions from Chat GPT)
% Original code from meeting in 07.2024
% Modified by: Jakub Nikonowicz since 08.2024 onwards
%
% Description: Configuration of SRS (Sounding Reference Signal) parameters
% as defined in 3GPP TS 38.211 Section 6.4.1.4.
%
%--------------------------------------------------------------------------
%
% Input parameters:
%   numSrsSymbols   - Number of OFDM symbols allocated per slot for SRS 
%                       transmission (1, 2, or 4)
%   symbolStartIdx  - Index of the starting OFDM symbol (0 to 13) within the 
%                       14-symbol OFDM slot
%   cSrs            - Bandwidth configuration index C_SRS (0...63), which 
%                       controls the allocated bandwidth for the SRS
%   combNumber      - Comb number (KTC), specifies frequency density in 
%                       subcarriers (2 or 4)
%--------------------------------------------------------------------------

    %% Define carrier and SRS configurations
    carrier = nrCarrierConfig;
    carrier.NSizeGrid = 273;           % Number of Resource Blocks; 273 = 100MHz
    carrier.SubcarrierSpacing = 30;    % Subcarrier spacing in kHz
    carrier.NSlot = 0;                 % Slot number for this example.
    
    %% Configure the SRS
    srs = nrSRSConfig;
    srs.NumSRSPorts = 1;               % Number of SRS antenna ports (1, 2, 4). 1 for UL Estimation for 1T4R typical UE.
    srs.SymbolStart = symbolStartIdx;  % Starting OFDM symbol within a slot.
    srs.NumSRSSymbols = numSrsSymbols; % Number of OFDM symbols allocated per slot (1, 2, 4).
    srs.ResourceType = 'periodic';     % Resource type ('periodic', 'semi-persistent', 'aperiodic'). Use 'aperiodic' to disable inter-slot frequency hopping.
    srs.SRSPeriod = [20 0];            % Periodicity and offset in slots. SRSPeriod(2) must be < SRSPeriod(1).
    srs.FrequencyStart = 0;            % Frequency position of the SRS in BWP in RBs.
    srs.NRRC = 0;                      % Additional offset from FreqStart specified in blocks of 4 PRBs (0...67).
    srs.CSRS = cSrs;                   % Bandwidth configuration C_SRS (0...63). Row index of bandwidth configuration Table 6.4.1.4.3-1: SRS bandwidth configuration.
    srs.BSRS = 0;                      % Bandwidth configuration B_SRS (0...3). Column index of bandwidth configuration table. Increasing the BSRS value decreases the SRS bandwidth.
    srs.BHop = 0;                      % Frequency hopping configuration (0...3). Set BHop < BSRS to enable frequency hopping.
    srs.KTC = combNumber;              % Comb number (2, 4), specifies frequency density in subcarriers.
    
    srs.SRSPositioning = 1;
    %srs
    
    %% Generate SRS indices and symbols
    [srsIndices, ~] = nrSRSIndices(carrier, srs);
    srsSymbols = nrSRS(carrier, srs);
    
    %% Print sizes of generated indices and symbols
    % disp('Size of srsIndices [(12 subcarrier / 2 KTC) x 272 PRB = 1632]:');
    % disp(size(srsIndices));
    % disp('Size and type of srsSymbols [compex double]:');
    % disp(size(srsSymbols));
    % disp(class(srsSymbols));
    
    %% Create and populate resource grid
    K = carrier.NSizeGrid * 12; % Total number of subcarriers
    L = 14; % Number of OFDM symbols per slot
    P = srs.NumSRSPorts; % Number of ports
    resourceGrid = complex(zeros(K, L, P));
    resourceGrid(srsIndices) = srsSymbols;
    
    %% Extract frequency coefficients 
    txFreqCoefficients = resourceGrid;  % 3276 active subcarriers
    % freqCoefficients is now a 3276 x L x P matrix
    
    %% OFDM Modulation Parameters
    ofdmInfo = nrOFDMInfo(carrier, 'CarrierFrequency', 3.5e9);
    txWaveForm = nrOFDMModulate(carrier, resourceGrid);
   
    %% Frequency domain analysis
    % % Define the sampling rate based on OFDM parameters
    % Fs = ofdmInfo.SampleRate;
    % N = length(txWaveform);
    % freqVector = (-N/2:N/2-1)*(Fs/N); % Frequency vector for FFT
    %
    % % Perform FFT and shift zero frequency component to the center
    % txWaveformFFT = fftshift(fft(txWaveform));
    % 
    % % Plot the amplitude of the subcarriers
    % figure;
    % plot(freqVector, abs(txWaveformFFT));
    % xlabel('Frequency (Hz)');
    % ylabel('Amplitude');
    % title('Amplitude of Subcarriers in Frequency Domain');
    % grid on;
    % 
    % % Plot the phase of the subcarriers
    % figure;
    % plot(freqVector, angle(txWaveformFFT));
    % xlabel('Frequency (Hz)');
    % ylabel('Phase (radians)');
    % title('Phase of Subcarriers in Frequency Domain');
    % grid on;

end
