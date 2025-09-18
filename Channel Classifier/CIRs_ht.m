function [h,index, delay_diff, pdp, label_output, H_time, H_delay, LoS_delay, max_Delay_new] = CIRs(tx_position, rx_position, scenario, rx_speed)
%clc; %clear all; %close all; 

%% Set the simulation parameters
% l.simpar.center_frequency              % Other way to set the simulation
% parameters
s = qd_simulation_parameters;            % New simulation parameters
s.center_frequency = 3.5e9;              % Carrier frequency, i.e., 3.80e9 for FR1 or 28e9 for FR2;
s.sample_density = 4;                    % 4 samples per half-wavelength
s.use_absolute_delays = 1;               % Include delay of the LOS path

%% Define Antenna 
a = qd_arrayant('omni');                 % Generate Omni antenna
a.center_frequency = s.center_frequency; 

%% Environment Laypout
l = qd_layout(s);                        % New Layout with the defined simulation parameters

[~, no_tx] = size(tx_position);  
% Transmitter
l.no_tx = no_tx;                        % BSs
l.tx_position = tx_position;            % Position of BSs at (x, y, height) in meters

for n=1 : no_tx 
    l.tx_array(1, n) = copy (a);                   % assign antenna "a" to BSs
end

% Receiver
% Track
t = qd_track('linear', 0, 0);           % 
t.initial_position = rx_position;       % Start position
t.set_speed(rx_speed)
%%%
l.no_rx = 1;                             % 1 Receiver                 
%l.rx_position = rx_position;             % Position of the receiver at (x, y, height) in meters  
l.rx_array = copy (a);
l.rx_track = t;

%Visualize the defined layout
%l.visualize;



%% Define channel model (i.e., scenario)
l.rx_track.scenario = { scenario };             % Scenario
p = l.init_builder;                             % Create channel builders
gen_parameters( p );                            % Generate small-scale fading
cn = get_channels( p );                         % Generate channel coefficients

%nFFT = 4096;                                   % FFT Size
BW = 3276*30e3;                                 % Bandwidth -> %30e3 for 100e6 (FR1) or 120e3 for 400e6 (FR2)

for n=1 : no_tx
    h(:, n) = cn(1, n).fr( BW, 3276 ); % Freq.-domain channel, 100/400 MHz bandwidth, 3276 carriers,
end  

% channel_time = struct([]);
% for i = 1 : size(cn,2)
%     channel_time(1,i).hd = [squeeze(cn(i).coeff) squeeze(cn(i).delay)]; % Combine h and d into a single array
% end

H_time = [];
H_delay = [];
for i = 1 : size(cn,2)
    H_time(:, i) = squeeze(cn(i).coeff);
    H_delay(:, i) = squeeze(cn(i).delay); % Combine h and d into a single array
end



% %Visualize channel response in time domain for BS = 1
% dd = squeeze((cn(1,2).delay));
% hh = squeeze((cn(1,2).coeff));
% figure;
% stem(dd,abs(hh), 'o');
% %hold on; % Keep the figure for additional plotting
% % Connect Points with Lines
% %plot(dd,20*log10(abs(hh)), '-r', 'LineWidth', 2); % Red line connecting points
% % Labels and Title
% xlabel('Delayes');
% ylabel('recievied power');
% %grid on; % Add grid for better visualization
% % Display the plot
% %hold off;
%% Verifications 
pdp = abs(ifft(h,[],1).'); %10*log10(abs(ifft(h,[],1).').^2);          % Power-delay profile
ind_pdp = 0:length(pdp)-1;                       % Calculate delays
delays  = ind_pdp/BW;
LoS_delay = sqrt(sum((l.tx_position - l.rx_position).^2))/299792458; %norm(l.tx_position-l.rx_position)/299792458; % True delay
%[maxVal, linearIdx] = max(A(:))
% fprintf('\n LOS_delay \n');
% disp(LoS_delay)
[maxVal, maxInd] = max(pdp, [],2);
[~,index] = max(pdp, [],1);
%fprintf('\n time stamp of max pow.: \n');
%max_Delay = maxInd/BW;
max_Delay_new = delays(maxInd);
% disp(max_Delay_new)
fprintf('\n Delay difference\n');
delay_diff = max_Delay_new-LoS_delay;
disp(delay_diff)
%% labelling the channels

% Threshold in nanoseconds
threshold = 10e-9;

% Logical check: Output 0 if < 10 ns, otherwise 1
label_output = delay_diff  >= threshold;

% Display the result
disp('Theoretical output based on delay difference:')
disp(label_output)



% %% Prepare Excel sheet
% xlsFilename = 'NLOS.xlsx';
% 
% % Check if file exists
% if isfile(xlsFilename)
%     % If file exists, count the number of existing sheets
%     [~, xlsSheets] = xlsfinfo(xlsFilename);
%     xlsSheetNumber = numel(xlsSheets) + 1; % New sheet numbered sequentially
% else
%     % If the file does not exist, we start from sheet 1
%     xlsSheetNumber = 1;
% end
% % Sheet name
% xlsSheetName = ['los_' num2str(xlsSheetNumber)];
% 
% % Adding headers
% xlsHeaders = {'Reference', 'LOS_diff', 'max_delay', 'los_delay', 'pdp'};
% writecell(xlsHeaders, xlsFilename, 'Sheet', xlsSheetName, 'Range', 'A1:J1');
% % We initialize the row variable to track the row number after headers
% xlsRow = 2;
% writematrix(delay_diff, xlsFilename, 'Sheet', xlsSheetName, 'Range', ['A' num2str(xlsRow)]);
% writematrix(max_Delay_new, xlsFilename, 'Sheet', xlsSheetName, 'Range', ['B' num2str(xlsRow)]);
% writematrix(LoS_delay, xlsFilename, 'Sheet', xlsSheetName, 'Range', ['C' num2str(xlsRow)]);
% writematrix(pdp, xlsFilename, 'Sheet', xlsSheetName, 'Range', ['D' num2str(xlsRow)]);

% figure(1)
% plot(delays(1:256), pdp(3,1:256)); hold on;
% xlabel('delay [s]');
% ylabel('power [dB]');
% title ( strvcat ( cn.name ) ); 
% % LOS_NLOS visualization

% %Plot the histogram (PDF)
% figure;
% histogram(delay_diff(:), 'Normalization', 'pdf'); % Normalize to PDF
% hold on;
% title('PDF and CDF of First Column');
% xlabel('Value');
% ylabel('Probability');
% grid on;
% title ( scenario );

% % Add Kernel Density Estimate for PDF
% [f, xi] = ksdensity(delay(:)); % Kernel density estimation
% plot(xi, f, 'r-', 'LineWidth', 2); % Overlay the smooth PDF

% % Compute and Plot CDF
% [cdfValues, cdfXi] = ecdf(delays(:)); % Empirical CDF
% plot(cdfXi, cdfValues, 'b-', 'LineWidth', 2); % Overlay the CDF
% 
% % Add legend
% % legend('Histogram (PDF)', 'CDF');
% % hold off;
% 
% % pdf_cdf
% 
% figure;
% 
% Plot the histogram (PDF approximation) on the left y-axis
% yyaxis left;
% histogram(delay_diff(:), 'Normalization', 'pdf');
% ylabel('Probability Density (PDF)');
% title ( scenario );
% xlabel('Delay Difference');
% grid on;
% 
% Compute the CDF
% [cdfValues, cdfXi] = ecdf(delay_diff(:));
% 
% Plot the CDF on the right y-axis
% yyaxis right;
% plot(cdfXi, cdfValues, 'r-', 'LineWidth', 2); % CDF curve
% ylabel('Cumulative Probability (CDF)');
% 
% Adjust legend
% legend('Histogram (PDF)', 'CDF', 'Location', 'best');
