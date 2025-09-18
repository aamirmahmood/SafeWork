clc;clear all;close all

d = load('rx_time_Data.mat');
%% 
rng(123)
% Initialize cell arrays for sequences and labels
data = d.data;

%% 
X = {};  % Input sequences
Y = {};  % One-hot labels

idx = 1;
for i = 1:8
    for j = 1:2
        for k = 1:50
            sample = data(i, j, k);
            % Transpose to [features × timesteps]
            X{idx, 1} = abs(sample.rx_time);  % Now it's [3 × 3276] per sample
            Y{idx, 1} = sample.labels';   % Column vector [3 × 1]
            idx = idx + 1;
        end
    end
end
% 
%% 
% Assume C is your 1×800 cell array
% labels = string(size(Y, 2));  % Preallocate output matrix
% 
% for i = 1:size(Y, 2)
%     row = Y{i}(:)';
%     if isequal(row, [false, false, false])
%         labels(i) = 'LOS';
%     elseif isequal(row, [true, true, true])
%         labels(i) = 'NLOS';
%     else
%         labels(i) = 'Semi';
%     end
% end
% 
% labels = categorical(labels);

XX = cell(size(X,1)*3,1);
counter = 1;
for i = 1:size(X,1)
    current_matrix = X{i};  % (3276x3)
    for j = 1:3
        XX{counter} = current_matrix(:, j); % Take each column separately
        counter = counter + 1;
    end
end

YMatrix = cell2mat(Y');
Y_temp = YMatrix(:);
YY = categorical(Y_temp);
% XMatrix = [];
% Y_onehot = [YMatrix == 0, YMatrix == 1];
% for i = 1:numel(X)
%     XMatrix = [XMatrix; X{i}];  % Stack the 3x3276 blocks vertically
% end
% X_cell = mat2cell(XMatrix, ones(size(XMatrix,1), 1), size(XMatrix, 2));
%% 
%% 
% % % Find minority class indices
% idxMinor = labels == "Semi-LOS/NLOS";
% XMinor = X(idxMinor)';
% YMinor = labels(idxMinor)';
% X = X(:);
% XMinor = XMinor(:);
% labels = labels(:);
% YMinor = YMinor(:);
% % Duplicate minority class
% X = [X; XMinor; XMinor]';
% labels = [labels; YMinor; YMinor];
% n = numel(labels);  % or size(XTrain, 1)
% % Generate random permutation of indices
% idx = randperm(n);
% % Shuffle both
% X = X(idx);
% labels = labels(idx);
%% 
%% 
% Assume: Y_onehot is [2400 × 2] logical one-hot matrix
% Define class names in order
% Y_onehot = double(Y_onehot);  % or single(Y_onehot)
% classNames = ["LOS", "NLOS"];
% YCat = onehotdecode(Y_onehot, classNames, 2); 

% Partition into 70% train, 30% temp
cv = cvpartition(YY, 'HoldOut', 0.3);
% Get training and test indices
idxTrain = training(cv);  % logical index vector
idxTest  = test(cv);      % logical index vector

XTrain = XX(idxTrain);
YTrain = YY(idxTrain);

XTemp = XX(idxTest);
YTemp = YY(idxTest);

% Split XTemp (30%) into 15% val, 15% test
cv2 = cvpartition(YTemp, 'HoldOut', 0.5);
idxTrain2 = training(cv2);  % logical index vector
idxTest2  = test(cv2);      % logical index vector

XVal = XTemp(idxTrain2);
YVal = YTemp(idxTrain2);

XTest = XTemp(idxTest2);
YTest = YTemp(idxTest2);

% Display sizes
disp("Training size: " + num2str(size(XTrain)));
disp("Validation size: " + num2str(size(XVal)));
disp("Test size: " + num2str(size(XTest)));
% 
% %% 
% for i = 1:numel(XTrain)
%     x = XTrain{i};
%     x = (x - mean(x)) ./ std(x);  % normalize
%     XTrain{i} = reshape(x, [], 1);  % ensure [3276 × 1]
% end
% 
% for i = 1:numel(XVal)
%     x = XVal{i};
%     x = (x - mean(x)) ./ std(x);
%     XVal{i} = reshape(x, [],1);
% end
% 
% for i = 1:numel(XTest)
%     x = XTest{i};
%     x = (x - mean(x)) ./ std(x);
%     XTest{i} = reshape(x, [],1);
% end
%% 

layers = [
    sequenceInputLayer(1, 'Normalization', 'zscore', 'MinLength', 3276)                              % 1 feature per time step
    %featureInputLayer()
    convolution1dLayer(16, 64, 'Stride', 2, 'Padding','same')
    batchNormalizationLayer
    reluLayer
    %dropoutLayer(0.4)

    % transposedConv1dLayer(8, 64, 'Stride', 2, 'Cropping', 'same')
    % batchNormalizationLayer
    % reluLayer

    convolution1dLayer(8, 32, 'Stride', 2, 'Padding','same')
    batchNormalizationLayer
    reluLayer
    %dropoutLayer(0.5)

    % transposedConv1dLayer(4, 32, 'Stride', 2, 'Cropping', 'same')
    % batchNormalizationLayer
    % reluLayer
    %dropoutLayer(0.4)

    globalAveragePooling1dLayer

    % fullyConnectedLayer(128)
    % reluLayer

    % 
    fullyConnectedLayer(64)
    reluLayer

    fullyConnectedLayer(32)
    reluLayer

    fullyConnectedLayer(2)
    softmaxLayer
    %classificationLayer
];

options = trainingOptions('adam', ...
    'InitialLearnRate', 1e-3, ...              % Starting learning rate
    'LearnRateSchedule', 'piecewise', ...      % Drop learning rate manually
    'LearnRateDropPeriod', 5, ...              % Drop every 5 epochs
    'LearnRateDropFactor', 0.6, ...            % Drop by half
    'MaxEpochs', 50, ...                       % Max number of epochs
    'MiniBatchSize', 64, ...                   % Batch size
    'Shuffle', 'every-epoch', ...              % Shuffle data each epoch
    'ValidationData', {XVal, YVal}, ...        % Validation set
    'ValidationFrequency', 10, ...             % Validate every 30 iterations
    'Verbose', true, ...                       % Print output to command line
    'VerboseFrequency', 30, ...                % Print every 30 iterations
    'ExecutionEnvironment', 'auto', ...
    "Metrics", ["accuracy", "auc"], ...         % Use GPU if available
    'GradientThreshold', 1, ...
    'ValidationPatience', 10, ...
    'L2Regularization', 1e-4, ...
    'Plots', 'training-progress')
% Train the network using trainnet 
[net,info] = trainnet(XTrain, YTrain, layers, "binary-crossentropy", options)
%% 
% %% 
% %YTest = renamecats(YTest, classNames);
% YPred = minibatchpredict(net,XTest);
% % Assume YPred is 120x3
% [~, classIdx] = max(YPred, [], 2);     % Get index of max in each row
% %YPred_hard = classIdx - 1
% % %% 
% % classNames = ["LOS", "NLOS", 'SEMI-LOS',]
% % predicted_labels = categorical(classNames(classIdx));
% % confusionchart(YTest, predicted_labels);
% %% 
% 
% % Initialize binary matrix
% YPred_bin = zeros(size(YPred));      % 120x3
% 
% % Set 1 at max index for each row
% for i = 1:size(YPred,1)
%     YPred_bin(i, classIdx(i)) = 1;
% end
% 
% % Step 1: Make sure both are column vectors
% YPred_label = YPred_bin(:);   % convert to 120×1
% YTest = YTest(:);               % ensure 120×1
% 
% % Step 2: Convert strings to categorical (match YTest categories)
% YPred_cat = categorical(YPred_label, categories(YTest));
% 
% % Step 3: Compute and visualize the confusion matrix
% confusionchart(YTest, YPred_cat);
% title('Confusion Matrix: YTest vs YPred');
% %% 
% classNames = ["LOS", "NLOS"];
% %YPred = onehotdecode(YPred,classNames,3);
% % Step 2: Absolute and normalized values
% YPred = minibatchpredict(net,XTest);
% [~, classIdx] = max(YPred, [], 2);
% binary_pred = categorical(classIdx)
% %binary_test = double(YTest == 'true');
% absVals = confusionmat(binary_test, classIdx);
% rowSums = sum(absVals, 2);
% normVals = absVals ./ rowSums;
%% 
classNames = ["LOS", "NLOS"];
YPred = minibatchpredict(net, XTest);
[~, classIdx] = max(YPred, [], 2);
% Convert predicted labels to binary
binary_classIdx = classIdx - 1;
% Ground truth
binary_test = double(YTest == 'true');
% Compute confusion matrix
absVals = confusionmat(binary_test, binary_classIdx);
% Normalize
rowSums = sum(absVals, 2);
normVals = absVals ./ rowSums;
%% 
% Step 3: Create confusion chart
fig = figure;
chart = confusionchart(binary_test, binary_classIdx, ...
    'Normalization', 'row-normalized')
        %'Title', 'Confusion Matrix with Absolute (Normalized) Values')
chart.FontSize = 15;
    drawnow
  % Ensure chart is fully rendered

% Step 4: Get axis limits to compute annotation positions
% Confusion matrix is 2x2 → grid coordinates: (1,1) to (2,2)
numClasses = size(absVals, 1);

% Step 5: Compute normalized positions for annotation placement
% Use the current axes position as reference
ax = chart.NormalizedValues;
pos = chart.Position;  % [left bottom width height]

for i = 1:numClasses
    for j = 1:numClasses
        % Create the label string
        txt = sprintf('%d (%.0f%%)', absVals(i, j), normVals(i, j)*100);

        % Compute annotation box coordinates
        x = pos(1) + (j - 0.1) * (pos(3) / numClasses);
        y = pos(2) + (numClasses - i + 0.1) * (pos(4) / numClasses); % invert y

        % Place annotation textbox
        annotation('textbox', [x-0.03, y-0.03, 0.06, 0.06], ...
            'String', txt, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'EdgeColor', 'none', ...
            'FontWeight', 'bold', ...
            'FontSize', 15);
    end
end
%% 
% Smooth with moving average (window size 10)
iterations = info.TrainingHistory{:,1};
loss = info.TrainingHistory{:,4}
accuracy = info.TrainingHistory{:,5}
auc = info.TrainingHistory{:,6};
smoothLoss = movmean(loss, 10);
smoothAcc = movmean(accuracy, 10);
smoothAUC = movmean(auc, 10);

figure;

% Plot 1: Loss
subplot(3,1,1);
plot(iterations, smoothLoss, 'b', 'LineWidth', 1.5); hold on;
plot(iterations, loss, ':k');
ylabel('Loss');
%title('Smoothed Training Loss');
legend('Smoothed', 'Raw');
grid on;

% Plot 2: Accuracy
subplot(3,1,2);
plot(iterations, smoothAcc, 'g', 'LineWidth', 1.5); hold on;
plot(iterations, accuracy, ':k');
ylabel('Accuracy (%)');
%title('Smoothed Training Accuracy');
legend('Smoothed', 'Raw');
grid on;

% Plot 3: AUC
subplot(3,1,3);
plot(iterations, smoothAUC, 'm', 'LineWidth', 1.5); hold on;
plot(iterations, auc, ':k');
xlabel('Iteration');
ylabel('AUC (%)');
%title('Smoothed Training AUC');
legend('Smoothed', 'Raw');
grid on;
% %% 
% layerSummary = analyzeNetwork(net);
% 
% % Extract Name and Description columns
% layerNames = layerSummary.Name;
% layerDesc  = layerSummary.Description;
% 
% % Combine into table
% T = table(layerNames, layerDesc, ...
%     'VariableNames', {'LayerName', 'Description'});
% 
% % Display table
% disp(T);
%%
% Convert to dlnetwork
%dlnet = dlnetwork(layers);
save('Conv1Modelf.mat', 'net');
