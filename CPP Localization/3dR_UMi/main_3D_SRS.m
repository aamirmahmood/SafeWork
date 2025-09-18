%Draft of a SRS-based scenario of CP in OFDM
clear all; clc; %#ok<CLALL>

%% Prepare Excel sheet
xlsFilename = 'simResultsSRS.xlsx';

% Check if file exists
if isfile(xlsFilename)
    % If file exists, count the number of existing sheets
    [~, xlsSheets] = xlsfinfo(xlsFilename);
    xlsSheetNumber = numel(xlsSheets) + 1; % New sheet numbered sequentially
else
    % If the file does not exist, we start from sheet 1
    xlsSheetNumber = 1;
end
% Sheet name
xlsSheetName = ['Sim_' num2str(xlsSheetNumber)];

% Adding headers
xlsHeaders = {'Reference', 'Estimated', 'Semi Results'};
writecell(xlsHeaders, xlsFilename, 'Sheet', xlsSheetName, 'Range', 'A1:J1');

% We initialize the row variable to track the row number after headers
xlsRow = 2;

%% Use "qd_builder.supported_scenarios" in command line to get all possible scenarios
% Other possible scenarios 
% 3GPP_38.901_InF_LOS
% 3GPP_38.901_InF_NLOS_SL
% 3GPP_38.901_InF_NLOS_SH
% 3GPP_38.901_InF_NLOS_DL
% 3GPP_38.901_InF_NLOS_DH

scenario = '3GPP_38.901_UMi_LOS'; % Urban microcells

%%RX movement (speed km/h)
rx_mov = 3;							

% Load data from CSV files
tx_data = readmatrix('antennas\bs_coord.csv', 'NumHeaderLines', 1); % BS coordinates
rx_data = readmatrix('antennas\new_input.csv', 'NumHeaderLines', 1); % Rx positions and BS indices
% Format: rx_x, rx_y, rx_z, BS_idx1, BS_idx2, BS_idx3

% Loop through receiver positions
for pos = 1:size(rx_data, 1)  % Iterate over each row in Rx data
    
    rx_position = rx_data(pos, 1:3)';   % Extract receiver position (x, y, z)
    bs_indices = rx_data(pos, 4:6);     % Extract indices of assigned BS (base stations)
    
    % Check if indices are within the valid range of BS data
    if all(bs_indices > 0 & bs_indices <= size(tx_data, 1))
        % Extract positions of the corresponding base stations
        tx_position(:,1) = tx_data(bs_indices(1), :)'; % BS-1 position
        tx_position(:,2) = tx_data(bs_indices(2), :)'; % BS-2 position
        tx_position(:,3) = tx_data(bs_indices(3), :)'; % BS-3 position
    else
        % Throw an error if indices are out of range
        error('Index out of range for BS data at position %d.', pos);
    end
    
    % Display the receiver and corresponding base station positions for this loop iteration
    disp(['Receiver position for pos ' num2str(pos) ':']);
    disp(rx_position');
    disp(['Base station positions for pos ' num2str(pos) ':']);
    disp(tx_position');

    % Dump reference position to the file
    writematrix(rx_position, xlsFilename, 'Sheet', xlsSheetName, 'Range', ['A' num2str(xlsRow)]);

    %% Generate a SRS block
    subcarrierSpacing = 30e3;
    numSrsSymbols = 8;  
    symbolStartIdx = 0;
    cSrs = 63;
    combNumber = 4; 
    
    [SRS_freq_tx, ~] = SRSSequence_5g(numSrsSymbols, symbolStartIdx, cSrs, combNumber);
    % SRS should be a 14-column vector of frequency coefficients. 273 RBs x 12 SC is 3276 subcarriers and for 30kHz SCS gives 100 MHz BW for each symbol.
    
    for symbolIdx = symbolStartIdx + 1 : symbolStartIdx + numSrsSymbols
    
        %% Get channel impulse response and create signal received at UE.
        [channel_freq_response, ~] = CIRs(tx_position, rx_position, scenario, rx_mov); % Channel impulse response modeled in QuaDriGa of size 3276 x N, where is N no. of gNbs's
        SRS_freq_rx = channel_freq_response .* SRS_freq_tx(:, symbolIdx);
    
        %% Checkpoint to show PRS at Tx and Rx
    %     figure(1)
    %         subplot(3,2,1); plot(abs(SRS_freq_tx(:, symbolIdx))); xlabel('subcarrier index [n]'); ylabel('amplitude SRS_{tx}');
    %         subplot(3,2,2); plot(angle(SRS_freq_tx(:, symbolIdx))); xlabel('subcarrier index [n]'); ylabel('angle SRS_{tx}');
    %         subplot(3,2,3); plot(abs(SRS_freq_rx)); xlabel('subcarrier index [n]'); ylabel('amplitude SRS_{rx}');
    %         subplot(3,2,4); plot(angle(SRS_freq_rx)); xlabel('subcarrier index [n]'); ylabel('angle SRS_{rx}');
    
         %% Prepare freq. block for processing 
         num_subcarriers = 3276;        % Actual no. of subcarriers
         speed_of_light = 299792458;    % Speed of light 
    
         % Excluding 0 freq. at Tx 
         % SRS_freq_nz_tx = []; 
         i=1;
         for n=1:num_subcarriers
            if SRS_freq_tx(n,symbolIdx) ~= 0
                SRS_freq_nz_tx(i) = SRS_freq_tx(n,symbolIdx); %#ok<*SAGROW> 
                i = i + 1;
            end
         end
    
     for j=1 : size(tx_position, 2)
    
         % Excluding 0 freq at Rx
         % SRS_freq_nz_rx = [];
         i=1;
         for n=1 : num_subcarriers
            if (SRS_freq_rx(n,j) ~= 0)
                SRS_freq_nz_rx(j,i) = SRS_freq_rx(n,j); 
                i = i + 1;
            end
         end
    
        % Phase shift by the initial offset.
        SRS_freq_nz_rx(j,:) = angle(SRS_freq_nz_rx(j,:)) - angle(SRS_freq_nz_tx);
        
        %% Differential phase measurements 
        % SRS_delta_phi_1 = [];
        % Multicase to keep the phase continuous 
        i=1;
        for n=2:size(SRS_freq_nz_rx, 2)
            if ((SRS_freq_nz_rx(j, n) - SRS_freq_nz_rx(j, n-1)) > pi)
                SRS_delta_phi_1(j, i) = abs((2*pi) - (SRS_freq_nz_rx(j, n) - SRS_freq_nz_rx(j, n-1)));
    
            elseif ((SRS_freq_nz_rx(j, n) - SRS_freq_nz_rx(j, n-1)) < -pi)
                SRS_delta_phi_1(j, i) = abs((2*pi) + (SRS_freq_nz_rx(j, n) - SRS_freq_nz_rx(j, n-1)));
    
            elseif ((SRS_freq_nz_rx(j, n) - SRS_freq_nz_rx(j, n-1)) < 0)
               SRS_delta_phi_1(j, i) = abs(-(SRS_freq_nz_rx(j, n) - SRS_freq_nz_rx(j, n-1)));
    
            else
                SRS_delta_phi_1(j, i) = abs((SRS_freq_nz_rx(j, n) - SRS_freq_nz_rx(j, n-1)));
            end
            i = i + 1;
        end
    
        delta_phi_1(j,:) = mean(SRS_delta_phi_1(j,:)./(2*pi));
        distance_1(j, symbolIdx) = delta_phi_1(j)/(1/(speed_of_light/(combNumber*subcarrierSpacing))); %/720000 -> 0.0024;
    
        %% Second Approach
        % SRS_delta_phi_m2 = [];
        i=1;
        for n=1:size(SRS_freq_nz_rx, 2)/2
            if ((SRS_freq_nz_rx(j, ((size(SRS_freq_nz_rx, 2)/2) + n)) - SRS_freq_nz_rx(j,n)) > pi)
                SRS_delta_phi_m2(j,i) = abs((2*pi) - (SRS_freq_nz_rx(j, (size(SRS_freq_nz_rx, 2)/2) + n) - SRS_freq_nz_rx(j,n)));
    
            elseif ((SRS_freq_nz_rx(j, ((size(SRS_freq_nz_rx, 2)/2) + n)) - SRS_freq_nz_rx(j,n)) < -pi)
                SRS_delta_phi_m2(j,i) = abs((2*pi) + (SRS_freq_nz_rx(j, (size(SRS_freq_nz_rx, 2)/2) + n) - SRS_freq_nz_rx(j,n)));
    
            elseif ((SRS_freq_nz_rx(j, ((size(SRS_freq_nz_rx, 2)/2) + n)) - SRS_freq_nz_rx(j,n)) < 0)
               SRS_delta_phi_m2(j,i) = abs(-((SRS_freq_nz_rx(j, (size(SRS_freq_nz_rx, 2)/2) + n)) - SRS_freq_nz_rx(j,n)));
    
            else
                SRS_delta_phi_m2(j,i) = abs((SRS_freq_nz_rx(j, (size(SRS_freq_nz_rx, 2)/2) + n)) - SRS_freq_nz_rx(j,n));
            end
            i = i + 1;
        end
    
        delta_phi_m2(j,:) = mean(SRS_delta_phi_m2(j,:)./(2*pi));
        del_N_m2(j) = floor((size(SRS_freq_nz_rx, 2)/2) * delta_phi_1(j) - delta_phi_m2(j));
        
        distance_m2(j, symbolIdx) = (delta_phi_m2(j) + del_N_m2(j))/(1/(speed_of_light/(num_subcarriers/2 * subcarrierSpacing))); %/196560000 -> 0.6552;
    
        %% Checkpoint to show results 
    %     figure(1) 
    %         subplot(3,2,5); plot(SRS_delta_phi_1(j,:)); xlabel('subcarrier index [n]'); ylabel('angle \Delta\phi_{1} [rad]');
    %         subplot(3,2,6); plot(SRS_delta_phi_m2(j,:)); xlabel('subcarrier index [n]'); ylabel('angle \Delta\phi_{m/2} [rad]');
    %     pause
    
     end
    end
    
    m_distance_1 = min(distance_1, [], 2);

    %% Trilateration
    p = tx_position;
    r = m_distance_1';
    [sol1, sol2] = trilateration(p, r);
    
    % Solutions
    fprintf('\n Solutions: \n');
    disp([sol1(2:4,1) sol2(2:4,1)])

    % Save estimated position to the file
    writematrix(sol2(2:4,1), xlsFilename, 'Sheet', xlsSheetName, 'Range', ['B' num2str(xlsRow)]);
    
        %% Checkpoint to show results 
    % Translate sphere to new location and plot as a surface.
%    figure(1)
%       [x1,y1,z1] = sphere;
%       [x2,y2,z2] = sphere;
%       [x3,y3,z3] = sphere;
%       % Scale to desire radius.
%       radius1 = m_distance_1(1,1);
%       radius2 = m_distance_1(2,1);
%       radius3 = m_distance_1(3,1);
%     
%       x1 = x1 * radius1;
%       y1 = y1 * radius1;
%       z1 = z1 * radius1;
%     
%       x2 = x2 * radius2;
%       y2 = y2 * radius2;
%       z2 = z2 * radius2;
%     
%       x3 = x3 * radius3;
%       y3 = y3 * radius3;
%       z3 = z3 * radius3;
%     
%     
%       surf(x1+tx_position(1,1), y1+tx_position(2,1), z1+tx_position(3,1), 'FaceAlpha',0.5); hold on
%       surf(x2+tx_position(1,2), y2+tx_position(2,2), z2+tx_position(3,2), 'FaceAlpha',0.5); hold on
%       surf(x3+tx_position(1,3), y3+tx_position(2,3), z3+tx_position(3,3), 'FaceAlpha',0.5);
%       % Label axes.
%       xlabel('X', 'FontSize', 5);
%       ylabel('Y', 'FontSize', 5);
%       zlabel('Z', 'FontSize', 5);
%       axis equal
%    pause
    
    %% Additional data
    semiResults = [];
    for c = 1:size(distance_1,2)
        rInt = distance_1(:,c)';
        [~, solInt] = trilateration(p, rInt);
        semiResults = [semiResults solInt(2:4,1)];
    end

    % Drop also semi-results to the file
    writematrix(semiResults, xlsFilename, 'Sheet', xlsSheetName, 'Range', ['C' num2str(xlsRow)]);

    %% Update row number
    xlsRow = xlsRow + 4;
end