%Draft of a PRS-based scenario of CP in OFDM
clear all; clc; %#ok<CLALL>

%% Prepare Excel sheet
xlsFilename = 'simResults.xlsx';

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
% 3GPP_38.901_UMi_LOS

%data = struct('rx_time', [], 'delay', [], 'labels', []);
data = struct('rx_time', [], 'delay_profile', [], 'labels', [], 'scenario', [], 'tx_coord', [], 'rx_coord', [], 'coord_est', [], ... 
    'dist_error', [], 'channel_coeff_time', [], 'true_delay', [],'est_delay', [], 'delay_difference', [], 'SRS_freq_rx', []);
Numtx = 3;
nFFT = 3276;
ML_flags = [];
net = load('Conv1Modelf.mat', 'net');
net = net.net;
for ii=1:2
    if ii == 1
        scenario = '3GPP_38.901_UMi_LOS';   % Channel model
    elseif ii == 2
        scenario = '3GPP_38.901_UMi_NLOS'; 
    else
        scenario = '3GPP_38.901_UMi_NLOS'; 
    end
%%%%%%%%%%%% tx BSs positions %%%%%%%%%%%%
    tx_position( ...
        : ,1) = [100; 100; 10];       % Position of BS-1 at (x; y; height) in meters
    tx_position( : ,2) = [150; 90; 10];       % Position of BS-2 at (x; y; height) in meters
    tx_position( : ,3) = [140; 150; 10]; 
    % tx_position( : ,4) = [0; 0; 18]; % Position of BS-3 at (x; y; height) in meters
    % tx_position( : ,5) = [80; 80; 18];
    % tx_position( : ,6) = [30; 120; 10];
    % tx_position( : ,7) = [100; 20; 10];
    % tx_position( : ,8) = [20; 100; 15];
    % tx_position( : ,9) = [100; 90; 10];
    % tx_position( : ,10) = [120; 80; 18]; % Position of BS-3 at (x; y; height) in meters
    
    rx_mov = [2.5; 5; 3];     %RX movement in (x; y; speed)
    
    rx_position_tot = [];
    est_tot = [];
    NumPos = 2;
    for pos = 0:NumPos-1
        rx_position = [120 + (rx_mov(1) *pos); 110 + (rx_mov(2) * pos); 1.5];
        rx_position_tot = [rx_position_tot, rx_position]; % Position of the receiver at (x; y; height) in meters
    
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
        delay = [];
        labels = [];
        rx_time = [];
        
        for symbolIdx = symbolStartIdx + 1 : symbolStartIdx + numSrsSymbols
        
            %% Get channel impulse response and create signal received at UE.
            [channel_freq_response, ~, delay_diff, pdp, label_output, H_time, H_delay, LoS_delay, max_Delay_new] = CIRs_ht(tx_position, rx_position, scenario, rx_mov(3));                  % Channel impulse response modeled in QuaDriGa of size 3276 x N, where is N no. of gNbs's
            %delay = [delay; H_delay];
            %labels = [labels; label_output];
            SRS_freq_rx = channel_freq_response .* SRS_freq_tx(:,symbolIdx);
            SRS_time_rx = ifft(SRS_freq_rx, nFFT, 1);
            
            bs_pred = [];
            for i=1:size(SRS_time_rx,2)
                YPred = minibatchpredict(net,{abs(SRS_time_rx(:,i)')'});
                [~, classIdx] = max(YPred, [], 2);
                binary_classIdx = classIdx - 1;
                %YPred = onehotdecode(YPred,classNames,2)
                bs_pred(end+1) = binary_classIdx;
            end
            ML_flags(end+1, :) = bs_pred;
            disp('Hey Jacob: My tiny ML predict that this SRS flags for three BSs is')
            disp(bs_pred)

            rx_time = cat(3,rx_time, SRS_time_rx); % concatenating the received signal for each position 
            delay = cat(3,delay, H_delay); % concatenating the delay spread for each position
            
            data(symbolIdx, ii, pos+1).rx_time = SRS_time_rx; %rx_time;
            data(symbolIdx, ii, pos+1).delay_profile = H_delay;
            data(symbolIdx, ii, pos+1).labels = label_output;
            data(symbolIdx, ii, pos+1).scenario = scenario;
            data(symbolIdx, ii, pos+1).rx_coord = rx_position;
            data(symbolIdx, ii, pos+1).channel_coeff_time = H_time;
            data(symbolIdx, ii, pos+1).true_delay = LoS_delay;
            data(symbolIdx, ii, pos+1).est_delay = max_Delay_new;
            data(symbolIdx, ii, pos+1).delay_difference = delay_diff;
            data(symbolIdx, ii, pos+1).SRS_freq_rx = SRS_freq_rx; %rx_freq;
            
            for j=1:Numtx
                data(symbolIdx, ii, pos+1).tx_coord(:,j) = tx_position(:,j);
            end

        end
    %end
    % 
    % 
%     % 
%     %         %% checkpoint check to see if received signal when is made in time domain is the same as frequency domain or not #### Verified ###
%     %         % NdelayTaps = 26;
%     %         % SRS_time_rx_check = ifft(SRS_freq_rx, NdelayTaps, 1);
%     %         % rx_time_check = zeros(size(SRS_time_rx_check)); % Preallocate the output matrix
%     %         % for i = 1:size(SRS_time_rx_check,2) % Iterate over each row
%     %         %     rx_time_check(:,i) = conv(SRS_time_rx_check(:,i), H_time(:,i), 'same'); % Convolution with 'same' size output
%     %         % end
%     %         %% Checkpoint to show PRS at Tx and Rx
%     %     %     figure(1)
%     %     %         subplot(3,2,1); plot(abs(SRS_freq_tx(:, symbolIdx))); xlabel('subcarrier index [n]'); ylabel('amplitude SRS_{tx}');
%     %     %         subplot(3,2,2); plot(angle(SRS_freq_tx(:, symbolIdx))); xlabel('subcarrier index [n]'); ylabel('angle SRS_{tx}');
%     %     %         subplot(3,2,3); plot(abs(SRS_freq_rx)); xlabel('subcarrier index [n]'); ylabel('amplitude SRS_{rx}');
%     %     %         subplot(3,2,4); plot(angle(SRS_freq_rx)); xlabel('subcarrier index [n]'); ylabel('angle SRS_{rx}');
%     % 
%     % %          %% Prepare freq. block for processing 
              num_subcarriers = 3276;        % Actual no. of subcarriers
              speed_of_light = 299792458;    % Speed of light 
    % 
    %          % Excluding 0 freq. at Tx 
    %          % SRS_freq_nz_tx = []; 
             i=1;
             for n=1:num_subcarriers
                if SRS_freq_tx(n,symbolIdx) ~= 0
                    SRS_freq_nz_tx(i) = SRS_freq_tx(n,symbolIdx); %#ok<*SAGROW> 
                    i = i + 1;
                end
             end
    % 
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
    % 
    %         % Phase shift by the initial offset.
            SRS_freq_nz_rx(j,:) = angle(SRS_freq_nz_rx(j,:)) - angle(SRS_freq_nz_tx);
    % 
    %         %% Differential phase measurements 
    %         % SRS_delta_phi_1 = [];
    %         % Multicase to keep the phase continuous 
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
    % 
            delta_phi_1(j,:) = mean(SRS_delta_phi_1(j,:)./(2*pi));
            distance_temp = delta_phi_1(j)/(1/(speed_of_light/(combNumber*subcarrierSpacing)));
            distance_1(j, symbolIdx) =  distance_temp;%/720000 -> 0.0024;

%     % 
%     %         %% Second Approach
%     %         % SRS_delta_phi_m2 = [];
%             i=1;
%             for n=1:size(SRS_freq_nz_rx, 2)/2
%                 if ((SRS_freq_nz_rx(j, ((size(SRS_freq_nz_rx, 2)/2) + n)) - SRS_freq_nz_rx(j,n)) > pi)
%                     SRS_delta_phi_m2(j,i) = abs((2*pi) - (SRS_freq_nz_rx(j, (size(SRS_freq_nz_rx, 2)/2) + n) - SRS_freq_nz_rx(j,n)));
% 
%                 elseif ((SRS_freq_nz_rx(j, ((size(SRS_freq_nz_rx, 2)/2) + n)) - SRS_freq_nz_rx(j,n)) < -pi)
%                     SRS_delta_phi_m2(j,i) = abs((2*pi) + (SRS_freq_nz_rx(j, (size(SRS_freq_nz_rx, 2)/2) + n) - SRS_freq_nz_rx(j,n)));
% 
%                 elseif ((SRS_freq_nz_rx(j, ((size(SRS_freq_nz_rx, 2)/2) + n)) - SRS_freq_nz_rx(j,n)) < 0)
%                    SRS_delta_phi_m2(j,i) = abs(-((SRS_freq_nz_rx(j, (size(SRS_freq_nz_rx, 2)/2) + n)) - SRS_freq_nz_rx(j,n)));
% 
%                 else
%                     SRS_delta_phi_m2(j,i) = abs((SRS_freq_nz_rx(j, (size(SRS_freq_nz_rx, 2)/2) + n)) - SRS_freq_nz_rx(j,n));
%                 end
%                 i = i + 1;
%             end
%     % 
%             delta_phi_m2(j,:) = mean(SRS_delta_phi_m2(j,:)./(2*pi));
%             del_N_m2(j) = floor((size(SRS_freq_nz_rx, 2)/2) * delta_phi_1(j) - delta_phi_m2(j));
%     % 
%             distance_m2(j, symbolIdx) = (delta_phi_m2(j) + del_N_m2(j))/(1/(speed_of_light/(num_subcarriers/2 * subcarrierSpacing))); %/196560000 -> 0.6552;
%     % % % 
%     % % %         %% Checkpoint to show results 
%     % % %         % figure(1) 
%     % % %         %     subplot(3,2,5); plot(SRS_delta_phi_1(j,:)); xlabel('subcarrier index [n]'); ylabel('angle \Delta\phi_{1} [rad]');
%     % % %         %     subplot(3,2,6); plot(SRS_delta_phi_m2(j,:)); xlabel('subcarrier index [n]'); ylabel('angle \Delta\phi_{m/2} [rad]');
%     % % %         % pause
%     % % %         % 
%          end
%          for j=1:symbolIdx
%             data(symbolIdx, ii, pos+1).coord_est = distance_1(:,symbolIdx);
%             data(symbolIdx, ii, pos+1).dist_error = sqrt(sum((rx_position - distance_1(:,symbolIdx)).^2));
          end
         
% 
%         %Writing the position estimations with respect to each symbols into
%         %the data structure
% 
%     % % % 
        m_distance_1 = min(distance_1, [], 2);
% 
%     % % % % 
%     % % % %     %% Trilateration
        p = tx_position;
        r = m_distance_1';
        [sol1, sol2] = trilateration(p, r);
        est_tot = [est_tot, sol2(2:4,1)];
        fprintf('\n est_values: \n');
        disp(est_tot)
        % Solutions
        fprintf('\n Solutions: \n');
        disp([sol1(2:4,1) sol2(2:4,1)])

        % Save estimated position to the file
        writematrix(sol2(2:4,1), xlsFilename, 'Sheet', xlsSheetName, 'Range', ['B' num2str(xlsRow)]);
%     % 
%     %         %% Checkpoint to show results 
%     %     % Translate sphere to new location and plot as a surface.
%     % %    figure(1)
%     % %       [x1,y1,z1] = sphere;
%     % %       [x2,y2,z2] = sphere;
%     % %       [x3,y3,z3] = sphere;
%     % %       % Scale to desire radius.
%     % %       radius1 = m_distance_1(1,1);
%     % %       radius2 = m_distance_1(2,1);
%     % %       radius3 = m_distance_1(3,1);
%     % %     
%     % %       x1 = x1 * radius1;
%     % %       y1 = y1 * radius1;
%     % %       z1 = z1 * radius1;
%     % %     
%     % %       x2 = x2 * radius2;
%     % %       y2 = y2 * radius2;
%     % %       z2 = z2 * radius2;
%     % %     
%     % %       x3 = x3 * radius3;
%     % %       y3 = y3 * radius3;
%     % %       z3 = z3 * radius3;
%     % % % %     
%     % % % %     
%     % % % %       surf(x1+tx_position(1,1), y1+tx_position(2,1), z1+tx_position(3,1), 'FaceAlpha',0.5); hold on
%     % % % %       surf(x2+tx_position(1,2), y2+tx_position(2,2), z2+tx_position(3,2), 'FaceAlpha',0.5); hold on
%     % % % %       surf(x3+tx_position(1,3), y3+tx_position(2,3), z3+tx_position(3,3), 'FaceAlpha',0.5);
%     % % % %       % Label axes.
%     % % % %       xlabel('X', 'FontSize', 5);
%     % % % %       ylabel('Y', 'FontSize', 5);
%     % % % %       zlabel('Z', 'FontSize', 5);
%     % % % %       axis equal
%     % % % %    pause
%     % % % 
%         %% Additional data
%         semiResults = [];
%         for c = 1:size(distance_1,2)
%             rInt = distance_1(:,c)';
%             [~, solInt] = trilateration(p, rInt);
%             semiResults = [semiResults solInt(2:4,1)];
%         end
% 
%         % Drop also semi-results to the file
%         writematrix(semiResults, xlsFilename, 'Sheet', xlsSheetName, 'Range', ['C' num2str(xlsRow)]);
% 
%         %% Update row number
%         xlsRow = xlsRow + 4;
%         end
% 
% % rx_time = num2cell(rx_time);
% % delay = num2cell(delay);
% % labels = num2cell(labels);
% % data(i+1, :) = {rx_time, delay, labels};
% %save('rx_time_LOS.mat', 'data');
% % %save('labels_LOS.mat', 'labels');
% % 
% %     data(i+1).Rx_Time = rx_time;
% %     data(i+1).Delay = delay;
% %     data(i+1).Labels = labels;
% %% Final estimation error
% 
%save('data.mat', 'data');
% % % % Compute the RMSE
rmse_2d = sqrt(mean((rx_position_tot(1:2,:) - est_tot(1:2,:)).^2, 'all'));
% % Display the result
disp('Root Mean Square Error (RMSE):');
disp(scenario);
disp(rmse_2d);
    end
end
 
% data = squeeze(reshape(data, [], NumPos));
%save('rx_time_Data.mat', 'data');

%% CHeck point to visualize LOS NLOS scenario label CDF 
% LOS_delay_scenario  = [];
% NLOS_delay_scenario = [];
% delay_diff = arrayfun(@(x) struct('delay_diff', x.delay_difference, 'scenarios', x.scenario), data(:));
% for i = 1:numel(delay_diff)
%     label_vec = delay_diff(i).delay_diff(:)';  % ensure 1x3 row vector
% 
%     if strcmp(hist_data(i).scenarios, "3GPP_38.901_UMi_LOS")
%         LOS_delay_scenario(end+1, :) = label_vec;     % append as a row
%     elseif strcmp(hist_data(i).scenarios, "3GPP_38.901_UMi_NLOS")
%         NLOS_delay_scenario(end+1, :) = label_vec;
%     end
% end
% flat_1 = abs(LOS_delay_scenario(:));
% flat_2 = abs(NLOS_delay_scenario(:));
% 
% % Combine and define logarithmic bins
% all_data = [flat_1; flat_2];
% min_val = floor(log10(min(all_data)));
% max_val = ceil(log10(max(all_data)));
% edges = logspace(min_val, max_val, 50);
% bin_centers = (edges(1:end-1) + edges(2:end)) / 2;
% 
% % Histogram and normalized CDFs
% counts_1 = histcounts(flat_1, edges);
% counts_2 = histcounts(flat_2, edges);
% 
% cdf_1 = cumsum(counts_1) / sum(counts_1);
% cdf_2 = cumsum(counts_2) / sum(counts_2);
% 
% % Plot CDFs
% figure;
% semilogx(bin_centers, cdf_1, 'r-', 'LineWidth', 2, 'DisplayName', 'LOS scenario');
% hold on;
% semilogx(bin_centers, cdf_2 , 'b--', 'LineWidth', 2, 'DisplayName', 'NLOS scenario');
% 
% xlabel('Value (log scale)');
% ylabel('Normalized CDF ');
% %title('Normalized CDFs with Logarithmic X-Axis');
% grid on;
% yticks([0 0.5 1 1+0.5 1+1])  % relabel y-axis
% yticklabels({'0','0.5','1','0.5','1'})
% legend('Location', 'best');