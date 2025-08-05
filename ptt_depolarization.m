clc; clear; close all;
addpath(genpath(pwd));
%%
parent_path = 'Dataset/Depolarization_Data';
bp_stim = 'cold';
scene = 'icu';
all_sbp = {};
all_dbp = {};
all_sid = {};
all_ptt = {};
sub_ids = 1:20;
data_path = fullfile(parent_path, upper(scene));
for enum_id = 1:length(sub_ids)
    sub_id = sub_ids(enum_id);
    for i_rec = 1:2
        data_fn = dir(fullfile(data_path, sprintf('%dP',sub_id),int2str(i_rec), '/*_win60s.mat'));
        load(fullfile(data_fn.folder, data_fn.name));
        load(fullfile(data_path, sprintf('%dP',sub_id),int2str(i_rec), '/PTT1_PTT2.mat'))
        min_length = min(length(SBP), length(PTT1));
        PTT1 = PTT1(:,1:min_length);
        SBP = SBP(1:min_length);
        DBP = DBP(1:min_length);
        all_sbp{end+1} = SBP;
        all_dbp{end+1} = DBP;
        all_ptt{end+1} = PTT1;
        all_sid{end+1} = (sub_id+10)*ones(1,min_length);
    end
end
%%
model_type = 'RF'; % 'SVR' or 'RF'
[all_Y_test_SBP, all_Y_test_DBP, all_Y_pred_SBP, all_Y_pred_DBP] = reg_bp(all_ptt, all_sbp, all_dbp, all_sid, model_type);
%%
figure;
plot(all_Y_test_SBP, 'k-', 'LineWidth', 1.5); hold on;
plot(all_Y_pred_SBP, 'r--', 'LineWidth', 1.5);
fill_between = fill([1:length(all_Y_test_SBP), fliplr(1:length(all_Y_test_SBP))], ...
    [all_Y_test_SBP(:)', fliplr(all_Y_pred_SBP(:)')], ...
    [1 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
legend('Ground Truth SBP', 'Predicted SBP', 'Prediction Error');
xlabel('Sample Index');
ylabel('SBP (mmHg)');
title('SBP Prediction with Error Region');
grid on; box on;
%%
figure;
plot(all_Y_test_DBP, 'k-', 'LineWidth', 1.5); hold on;
plot(all_Y_pred_DBP, 'r--', 'LineWidth', 1.5);
fill_between = fill([1:length(all_Y_test_DBP), fliplr(1:length(all_Y_test_DBP))], ...
    [all_Y_test_DBP(:)', fliplr(all_Y_pred_DBP(:)')], ...
    [1 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
legend('Ground Truth DBP', 'Predicted DBP', 'Prediction Error');
xlabel('Sample Index');
ylabel('DBP (mmHg)');
title('DBP Prediction with Error Region');
grid on; box on;
%%
% Extract subject list
subjects = unique(cell2mat(all_sid));
n_subjects = length(subjects);

% Initialize cursor
start_idx = 1;

% Determine subplot layout
n_rows = ceil(sqrt(n_subjects));
n_cols = ceil(n_subjects / n_rows);

figure('Name','SBP Prediction per Subject','Position',[100 100 1400 800]);

for i = 1:n_subjects
    sid = subjects(i);

    % Collect true subject data length
    this_len = 0;
    for j = 1:length(all_sid)
        if all_sid{j}(1) == sid
            this_len = this_len + size(all_ptt{j}, 2);  % number of timepoints
        end
    end

    idx_range = start_idx : (start_idx + this_len - 1);
    
    % Plot
    subplot(n_rows, n_cols, i);
    plot(all_Y_test_SBP(idx_range), 'k-', 'LineWidth', 1.2); hold on;
    plot(all_Y_pred_SBP(idx_range), 'r--', 'LineWidth', 1.2);
    title(sprintf('Subject %d', sid));
    xlabel('Sample Index'); ylabel('SBP');
    xlim([1 length(idx_range)]);
    grid on; box on;

    if i == 1
        legend('Ground Truth', 'Prediction', 'Location', 'best');
    end

    % Move to next subject range
    start_idx = start_idx + this_len;
end
%%
% Extract subject list
subjects = unique(cell2mat(all_sid));
n_subjects = length(subjects);

% Initialize cursor
start_idx = 1;

% Determine subplot layout
n_rows = ceil(sqrt(n_subjects));
n_cols = ceil(n_subjects / n_rows);

figure('Name','DBP Prediction per Subject','Position',[100 100 1400 800]);

for i = 1:n_subjects
    sid = subjects(i);

    % Collect true subject data length
    this_len = 0;
    for j = 1:length(all_sid)
        if all_sid{j}(1) == sid
            this_len = this_len + size(all_ptt{j}, 2);  % number of timepoints
        end
    end

    idx_range = start_idx : (start_idx + this_len - 1);
    
    % Plot
    subplot(n_rows, n_cols, i);
    plot(all_Y_test_DBP(idx_range), 'k-', 'LineWidth', 1.2); hold on;
    plot(all_Y_pred_DBP(idx_range), 'r--', 'LineWidth', 1.2);
    title(sprintf('Subject %d', sid));
    xlabel('Sample Index'); ylabel('DBP');
    xlim([1 length(idx_range)]);
    grid on; box on;

    if i == 1
        legend('Ground Truth', 'Prediction', 'Location', 'best');
    end

    % Move to next subject range
    start_idx = start_idx + this_len;
end
%%
% Setup
subjects = unique(cell2mat(all_sid));
% n_subjects = length(subjects);
n_subjects = 5;
start_idx = 1;

figure('Name','SBP Prediction and PTT1 per Subject','Position',[100 100 1400 2*n_subjects*100]);

for i = 1:n_subjects
    sid = subjects(i);

    % Collect data for this subject
    this_len = 0;
    this_ptt = [];  % 3 × T
    for j = 1:length(all_sid)
        if all_sid{j}(1) == sid
            this_len = this_len + size(all_ptt{j}, 2);
            this_ptt = [this_ptt, all_ptt{j}];  % concatenate along time
        end
    end
    idx_range = start_idx : (start_idx + this_len - 1);

    % Row 1: SBP Prediction
    subplot(n_subjects, 3, 3*(i-1)+1);
    plot(all_Y_test_SBP(idx_range), 'k-', 'LineWidth', 1.2); hold on;
    plot(all_Y_pred_SBP(idx_range), 'r--', 'LineWidth', 1.2);
    title(sprintf('Subject %d - SBP Prediction', sid));
    xlabel('Sample Index'); ylabel('SBP (mmHg)');
    xlim([1 length(idx_range)]);
    grid on; box on;
    if i == 1
        legend('Ground Truth', 'Prediction', 'Location', 'best');
    end

    % Row 2: DBP Prediction
    subplot(n_subjects, 3, 3*(i-1)+2);
    plot(all_Y_test_DBP(idx_range), 'k-', 'LineWidth', 1.2); hold on;
    plot(all_Y_pred_DBP(idx_range), 'r--', 'LineWidth', 1.2);
    title(sprintf('Subject %d - DBP Prediction', sid));
    xlabel('Sample Index'); ylabel('DBP (mmHg)');
    xlim([1 length(idx_range)]);
    grid on; box on;
    if i == 2
        legend('Ground Truth', 'Prediction', 'Location', 'best');
    end

    % Row 3: PTT1 Channels
    subplot(n_subjects, 3, 3*(i-1)+3);
    time_vec = 1:size(this_ptt, 2);
    plot(time_vec, this_ptt(1,:), 'b', ...
         time_vec, this_ptt(2,:), 'g', ...
         time_vec, this_ptt(3,:), 'm', 'LineWidth', 1);
    title(sprintf('Subject %d - PTT1 (3 Channels)', sid));
    xlabel('Sample Index'); ylabel('PTT Value');
    if i == 3
        legend('Ch1','Ch2','Ch3','Location','best');
    end
    
    xlim([1 size(this_ptt,2)]);
    grid on; box on;

    % Update index
    start_idx = start_idx + this_len;
end
%%
% Setup
subjects = unique(cell2mat(all_sid));
% n_subjects = length(subjects);
n_subjects = 5;
start_idx = 1;

figure('Name','DBP Prediction and PTT1 per Subject','Position',[100 100 2000 3*n_subjects*100]);

for i = 1:n_subjects
    sid = subjects(i);

    % Collect data for this subject
    this_len = 0;
    this_ptt = [];  % 3 × T
    for j = 1:length(all_sid)
        if all_sid{j}(1) == sid
            this_len = this_len + size(all_ptt{j}, 2);
            this_ptt = [this_ptt, all_ptt{j}];  % concatenate along time
        end
    end
    idx_range = start_idx : (start_idx + this_len - 1);

    % Row 1: DBP Prediction
    subplot(n_subjects, 2, 2*(i-1)+1);
    plot(all_Y_test_DBP(idx_range), 'k-', 'LineWidth', 1.2); hold on;
    plot(all_Y_pred_DBP(idx_range), 'r--', 'LineWidth', 1.2);
    title(sprintf('Subject %d - DBP Prediction', sid));
    xlabel('Sample Index'); ylabel('DBP (mmHg)');
    xlim([1 length(idx_range)]);
    grid on; box on;
    if i == 1
        legend('Ground Truth', 'Prediction', 'Location', 'best');
    end

    % Row 2: PTT1 Channels
    subplot(n_subjects, 2, 2*(i-1)+2);
    time_vec = 1:size(this_ptt, 2);
    plot(time_vec, this_ptt(1,:), 'b', ...
         time_vec, this_ptt(2,:), 'g', ...
         time_vec, this_ptt(3,:), 'm', 'LineWidth', 1);
    title(sprintf('Subject %d - PTT1 (3 Channels)', sid));
    xlabel('Sample Index'); ylabel('PTT Value');
    legend('Ch1','Ch2','Ch3','Location','best');
    xlim([1 size(this_ptt,2)]);
    grid on; box on;

    % Update index
    start_idx = start_idx + this_len;
end
%%
% Config
model_type = 'SVR'; % 'SVR' or 'RF'

function [all_Y_test_SBP, all_Y_test_DBP, all_Y_pred_SBP, all_Y_pred_DBP] = reg_bp(all_ptt, all_sbp, all_dbp, all_sid, model_type)
% Config
rng(42);  % for reproducibility

% Extract subject list
subjects = unique(cell2mat(all_sid));
n_subjects = length(subjects);

% Initialize metrics
all_metrics = struct('SBP_MAE', [], 'SBP_RMSE', [], 'SBP_R', [], ...
                     'DBP_MAE', [], 'DBP_RMSE', [], 'DBP_R', []);
all_Y_test_SBP = []; all_Y_pred_SBP = [];
all_Y_test_DBP = []; all_Y_pred_DBP = [];

for i = 1:n_subjects
    test_sid = subjects(i);

    % Init train/test sets
    X_train = []; Y_train_SBP = []; Y_train_DBP = [];
    X_test  = []; Y_test_SBP  = []; Y_test_DBP  = [];

    for j = 1:length(all_sid)
        sid = all_sid{j}(1);
        PTT = all_ptt{j};       % 3 × T
        SBP = all_sbp{j};       % T × 1
        DBP = all_dbp{j};

        % Transpose PTT to T × 3
        X = PTT';  % T × 3
        y_sbp = SBP(:);
        y_dbp = DBP(:);

        if sid == test_sid
            X_test = [X_test; X];
            Y_test_SBP = [Y_test_SBP; y_sbp];
            Y_test_DBP = [Y_test_DBP; y_dbp];
        else
            X_train = [X_train; X];
            Y_train_SBP = [Y_train_SBP; y_sbp];
            Y_train_DBP = [Y_train_DBP; y_dbp];
        end
    end

    % Train models
    if strcmpi(model_type, 'SVR')
        mdl_SBP = fitrsvm(X_train, Y_train_SBP, 'KernelFunction', 'rbf');
        mdl_DBP = fitrsvm(X_train, Y_train_DBP, 'KernelFunction', 'rbf');
    elseif strcmpi(model_type, 'RF')
        mdl_SBP = fitrensemble(X_train, Y_train_SBP, 'Method', 'Bag');
        mdl_DBP = fitrensemble(X_train, Y_train_DBP, 'Method', 'Bag');
    else
        error('Unknown model type.');
    end

    % Predict
    Y_pred_SBP = predict(mdl_SBP, X_test);
    Y_pred_DBP = predict(mdl_DBP, X_test);

    % Append to overall record
    all_Y_test_SBP = [all_Y_test_SBP; Y_test_SBP];
    all_Y_pred_SBP = [all_Y_pred_SBP; Y_pred_SBP];
    all_Y_test_DBP = [all_Y_test_DBP; Y_test_DBP];
    all_Y_pred_DBP = [all_Y_pred_DBP; Y_pred_DBP];

    % Per-subject metrics
    sbp_mae = mean(abs(Y_pred_SBP - Y_test_SBP));
    sbp_rmse = sqrt(mean((Y_pred_SBP - Y_test_SBP).^2));
    sbp_r = corr(Y_pred_SBP, Y_test_SBP);

    dbp_mae = mean(abs(Y_pred_DBP - Y_test_DBP));
    dbp_rmse = sqrt(mean((Y_pred_DBP - Y_test_DBP).^2));
    dbp_r = corr(Y_pred_DBP, Y_test_DBP);

    all_metrics.SBP_MAE(end+1) = sbp_mae;
    all_metrics.SBP_RMSE(end+1) = sbp_rmse;
    all_metrics.SBP_R(end+1) = sbp_r;

    all_metrics.DBP_MAE(end+1) = dbp_mae;
    all_metrics.DBP_RMSE(end+1) = dbp_rmse;
    all_metrics.DBP_R(end+1) = dbp_r;

    fprintf('Subject %2d | SBP MAE=%.2f RMSE=%.2f R=%.2f | DBP MAE=%.2f RMSE=%.2f R=%.2f\n', ...
        test_sid, sbp_mae, sbp_rmse, sbp_r, dbp_mae, dbp_rmse, dbp_r);
end

% Overall pooled metrics
overall_mae_sbp = mean(abs(all_Y_pred_SBP - all_Y_test_SBP));
overall_rmse_sbp = sqrt(mean((all_Y_pred_SBP - all_Y_test_SBP).^2));

overall_mae_dbp = mean(abs(all_Y_pred_DBP - all_Y_test_DBP));
overall_rmse_dbp = sqrt(mean((all_Y_pred_DBP - all_Y_test_DBP).^2));

% Final report
fprintf('\n=== Leave-One-Subject-Out Mean Results ===\n');
fprintf('SBP: MAE = %.2f ± %.2f | RMSE = %.2f ± %.2f | R = %.2f ± %.2f\n', ...
    mean(all_metrics.SBP_MAE), std(all_metrics.SBP_MAE), ...
    mean(all_metrics.SBP_RMSE), std(all_metrics.SBP_RMSE), ...
    mean(all_metrics.SBP_R), std(all_metrics.SBP_R));

fprintf('DBP: MAE = %.2f ± %.2f | RMSE = %.2f ± %.2f | R = %.2f ± %.2f\n', ...
    mean(all_metrics.DBP_MAE), std(all_metrics.DBP_MAE), ...
    mean(all_metrics.DBP_RMSE), std(all_metrics.DBP_RMSE), ...
    mean(all_metrics.DBP_R), std(all_metrics.DBP_R));

fprintf('\n=== Overall Pooled Prediction Statistics ===\n');
fprintf('SBP: MAE = %.2f | RMSE = %.2f\n', overall_mae_sbp, overall_rmse_sbp);
fprintf('DBP: MAE = %.2f | RMSE = %.2f\n', overall_mae_dbp, overall_rmse_dbp);
end
%%
