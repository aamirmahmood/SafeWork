%Draft of a PRS-based scenario of CP in OFDM
clear all; clc; %#ok<CLALL>

%% Prepare Excel sheet
xlsFilename = 'simResultsPRS.xlsx';

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
rx_data = readmatrix('antennas\new_input_test.csv', 'NumHeaderLines', 1); % Rx positions and BS indices
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

    for trial_cnt = 1:8
     
        %% Generate a PRS block
        fc = 3.5e9;  %28e9 or 3.8e9
        fra  = 30e3; %120e3 or 30e3
        Ts = 1/(fra*4096);
        NPRB = 273;
        l1   = 1;
        ns   = 1;
        beta_PRS = 1;
        lstart = 1;
        BSID = 1;
        PRS_freq_tx = PRSSequence_5g( NPRB,l1,BSID+100,ns,beta_PRS,lstart );
        % PRS should be a single-column vector of frequency coefficients for 3276 subcarriers. 100MHz -> 273 RB x 12 = 3276 SC;
    
        %% Create reference signal at gNB.
        % Reference signal at UE. Loop over symbols, aggregated in time.
        [nFFT,number_of_symbols] = size(PRS_freq_tx);
        PRS_time_tx = [];
        for k_symbol = 1:number_of_symbols
            PRS_time_tx = [PRS_time_tx sqrt(nFFT)*ifft(PRS_freq_tx(:,k_symbol),nFFT)']; %#ok<AGROW> % PRS_freq_tx should be a matrix with ifft done in columns
        end
    
        %% Get channel impulse response and create signal received at UE.
        [channel_freq_response, ~] = CIRs(tx_position, rx_position, scenario, rx_mov);	% Channel impulse response modeled in QuaDriGa of size 3276xN; N - no. of gNbs's
        channel_freq_response = padarray(channel_freq_response,4096-3276,0,'post');		% Channel impulse response should be a colon vector of size 4096x1  
        channel_freq_response = repmat(channel_freq_response,1,number_of_symbols);      % Debugging with channel_response = ones(4096,1);
        PRS_freq_rx = channel_freq_response .* PRS_freq_tx;
    
        %% Checkpoint to show PRS at Tx and Rx
%         if pos == 1 && trial_cnt == 1 
%             figure(1);
%             
%             % Plot the first subplot
%             subplot(2,1,1);
%             plot(angle(PRS_freq_tx(:,1)));
%             xlim([1 3276]);
%             xlabel('Subcarrier index'); 
%             ylabel('Phase [rad]');
%             legend('Phase spectrum at Tx');
%             
%             % Get handle for the first subplot
%             ax1 = gca;
%             
%             % Plot the second subplot
%             subplot(2,1,2);
%             plot(angle(PRS_freq_rx(:,1))); 
%             xlim([1 3276]);
%             xlabel('Subcarrier index'); 
%             ylabel('Phase [rad]');
%             legend('Phase spectrum at Rx');
%             
%             % Get handle for the second subplot
%             ax2 = gca;
%             
%             % Apply settings to both subplots
%             for ax = [ax1, ax2]
%                 set(ax, 'FontSize', 12);           % Set font size for axes
%                 set(ax, 'GridAlpha', 0.3);         % Grid transparency
%                 grid(ax, 'on');                    % Enable grid
%                 box(ax, 'on');                     % Enable box around the plot
%                 set(ax, 'Color', 'none');          % Set axes background to transparent
%             end
%             
%             % Set figure background to transparent
%             set(gcf, 'Color', 'none');
%             
%             % Export the figure to a PNG file
%             outputFile = 'fig_cpp_2.png';          % Output file name
%             print(outputFile, '-dpng', '-r600');   % Export to PNG with 600 DPI
%         end
    
         %% Prepare freq. block for processing 
         num_subcarriers = 3276;        % Actual no. of subcarriers
         speed_of_light = 299792458;    % Speed of light 
    
         % Excluding 0 freq. at Tx 
         PRS_freq_nz_tx = []; 
         i=1;
         for n=1:num_subcarriers
            if PRS_freq_tx(n) ~= 0
                PRS_freq_nz_tx(i) = PRS_freq_tx(n); %#ok<*SAGROW> 
                i = i + 1;
            end
         end
    
     for j=1 : size(tx_position, 2)
    
         % Excluding 0 freq at Rx
         %PRS_freq_nz_rx = [];
         i=1;
         for n=1 : num_subcarriers
            if (PRS_freq_rx(n,j) ~= 0)
                PRS_freq_nz_rx(j,i) = PRS_freq_rx(n,j); 
                i = i + 1;
            end
         end
    
        % Phase shift by the initial offset.
        PRS_freq_nz_rx_rect(j,:) = angle(PRS_freq_nz_rx(j,:)) - angle(PRS_freq_nz_tx);
        
        %% Checkpoint to show results 
%         if pos == 1 && trial_cnt == 1 
%             figure(2);
%             
%             % Plot the first subplot
%             subplot(2,1,1);
%             plot(angle(PRS_freq_rx(:,1)));
%             xlim([1 3276]);
%             xlabel('Subcarrier index'); 
%             ylabel('Phase [rad]');
%             legend('Phase spectrum at Rx');
%             
%             % Get handle for the first subplot
%             ax1 = gca;
%             
%             % Plot the second subplot
%             subplot(2,1,2);
%             plot(PRS_freq_nz_rx_rect(1,:)); 
%             xlabel('Subcarrier index'); 
%             ylabel('Phase [rad]');
%             legend(['Phase spectrum at Rx' newline 'offset-corrected' newline 'non-zero SCs']);
%             
%             % Get handle for the second subplot
%             ax2 = gca;
%             
%             % Apply settings to both subplots
%             for ax = [ax1, ax2]
%                 set(ax, 'FontSize', 12);           % Set font size for axes
%                 set(ax, 'GridAlpha', 0.3);         % Grid transparency
%                 grid(ax, 'on');                    % Enable grid
%                 box(ax, 'on');                     % Enable box around the plot
%                 set(ax, 'Color', 'none');          % Set axes background to transparent
%             end
%             
%             % Set figure background to transparent
%             set(gcf, 'Color', 'none');
%             
%             % Export the figure to a PNG file
%             outputFile = 'fig_cpp_3.png';          % Output file name
%             print(outputFile, '-dpng', '-r600');   % Export to PNG with 600 DPI
%         end        
        
        %% Differential phase measurements 
        %PRS_delta_phi_1 = [];
        % Multicase to keep the phase continuous 
        i=1;
        for n=2:size(PRS_freq_nz_rx_rect, 2)
            if ((PRS_freq_nz_rx_rect(j, n) - PRS_freq_nz_rx_rect(j, n-1)) > pi)
                PRS_delta_phi_1(j, i) = abs((2*pi) - (PRS_freq_nz_rx_rect(j, n) - PRS_freq_nz_rx_rect(j, n-1)));
    
            elseif ((PRS_freq_nz_rx_rect(j, n) - PRS_freq_nz_rx_rect(j, n-1)) < -pi)
                PRS_delta_phi_1(j, i) = abs((2*pi) + (PRS_freq_nz_rx_rect(j, n) - PRS_freq_nz_rx_rect(j, n-1)));
    
            elseif ((PRS_freq_nz_rx_rect(j, n) - PRS_freq_nz_rx_rect(j, n-1)) < 0)
               PRS_delta_phi_1(j, i) = abs(-(PRS_freq_nz_rx_rect(j, n) - PRS_freq_nz_rx_rect(j, n-1)));
    
            else
                PRS_delta_phi_1(j, i) = abs((PRS_freq_nz_rx_rect(j, n) - PRS_freq_nz_rx_rect(j, n-1)));
            end
            i = i + 1;
        end
    
        delta_phi_1(j,:) = mean(PRS_delta_phi_1(j,:)./(2*pi));
        distance_1(j, trial_cnt) = delta_phi_1(j)/(1/(speed_of_light/180000)); %0.0024; 720000 for 120e3

        temp_PRS_freq_nz_rx_rect = PRS_freq_nz_rx_rect(:,2:end) - PRS_freq_nz_rx_rect(:,1:end-1);
        
        %% Checkpoint to show results 
%         if pos == 1 && trial_cnt == 1 
%             figure(3);
%             
%             % Plot the first subplot
%             subplot(2,1,1);
%             plot(temp_PRS_freq_nz_rx_rect(1,:)); 
%             xlabel('Subcarrier index'); 
%             ylabel('\Delta\phi_{1} [rad]');
%             legend('NZ OC phase spectrum at Rx');
%             
%             % Get handle for the first subplot
%             ax1 = gca;
%             
%             % Plot the second subplot
%             subplot(2,1,2);
%             plot(PRS_delta_phi_1(1,:)); 
%             xlabel('Subcarrier index'); 
%             ylabel('\Delta\phi_{1} [rad]');
%             legend(['NZ OC phase spectrum at Rx' newline 'with phase progression']);
%             
%             % Get handle for the second subplot
%             ax2 = gca;
%             
%             % Apply settings to both subplots
%             for ax = [ax1, ax2]
%                 set(ax, 'FontSize', 12);           % Set font size for axes
%                 set(ax, 'GridAlpha', 0.3);         % Grid transparency
%                 grid(ax, 'on');                    % Enable grid
%                 box(ax, 'on');                     % Enable box around the plot
%                 set(ax, 'Color', 'none');          % Set axes background to transparent
%             end
%             
%             % Set figure background to transparent
%             set(gcf, 'Color', 'none');
%             
%             % Export the figure to a PNG file
%             outputFile = 'fig_cpp_4.png';          % Output file name
%             print(outputFile, '-dpng', '-r600');   % Export to PNG with 600 DPI
% 
%         end
    
       %% Second Approach
       % delta_phi_m2 = []; %#ok<NASGU> 
       % PRS_delta_phi_m2 = [];
        i=1;
        for n=1:size(PRS_freq_nz_rx_rect, 2)/2
            if ((PRS_freq_nz_rx_rect(j, ((size(PRS_freq_nz_rx_rect, 2)/2) + n)) - PRS_freq_nz_rx_rect(j,n)) > pi)
                PRS_delta_phi_m2(j,i) = abs((2*pi) - (PRS_freq_nz_rx_rect(j, (size(PRS_freq_nz_rx_rect, 2)/2) + n) - PRS_freq_nz_rx_rect(j,n)));
    
            elseif ((PRS_freq_nz_rx_rect(j, ((size(PRS_freq_nz_rx_rect, 2)/2) + n)) - PRS_freq_nz_rx_rect(j,n)) < -pi)
                PRS_delta_phi_m2(j,i) = abs((2*pi) + (PRS_freq_nz_rx_rect(j, (size(PRS_freq_nz_rx_rect, 2)/2) + n) - PRS_freq_nz_rx_rect(j,n)));
    
            elseif ((PRS_freq_nz_rx_rect(j, ((size(PRS_freq_nz_rx_rect, 2)/2) + n)) - PRS_freq_nz_rx_rect(j,n)) < 0)
               PRS_delta_phi_m2(j,i) = abs(-((PRS_freq_nz_rx_rect(j, (size(PRS_freq_nz_rx_rect, 2)/2) + n)) - PRS_freq_nz_rx_rect(j,n)));
    
            else
                PRS_delta_phi_m2(j,i) = abs((PRS_freq_nz_rx_rect(j, (size(PRS_freq_nz_rx_rect, 2)/2) + n)) - PRS_freq_nz_rx_rect(j,n));
            end
            i = i + 1;
        end
    
        delta_phi_m2(j,:) = mean(PRS_delta_phi_m2(j,:)./(2*pi));
        del_N_m2(j) = floor((size(PRS_freq_nz_rx_rect, 2)/2) * delta_phi_1(j) - delta_phi_m2(j));
        
        distance_m2(j, trial_cnt) = (delta_phi_m2(j) + del_N_m2(j))/(1/(speed_of_light/49140000)); %0.6552; 196560000 for 120e3
    
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
   