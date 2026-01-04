%==========================================================================
% SPAD Detection Simulation - Publication-Ready Visualization
%==========================================================================
% Description:
%   Generates high-quality timeline visualization of SPAD (Single Photon
%   Avalanche Diode) detection events across multiple counting periods
%   with varying photon rates. Demonstrates dead time effects and event
%   losses. Suitable for manuscript figures.
%
% Visualization Features:
%   - Color-coded periods and dead time regions
%   - Detected events (filled circles) and missed events (crosses)
%   - Period boundaries and timeline markers
%   - Vector output for publication quality
%
% Author: Akatsuki Sky
% Date: 2026-1-1
%==========================================================================

clear; clc;
rng(42);

%----------------------------------------------------------------------
% Simulation Parameters
%----------------------------------------------------------------------
Ts = 180;          % Period length in ns
tau = 45;          % Dead time in ns

mu_k_vec = [1, 10, 100]; % Expected event count per period (varying load)

numPeriods = numel(mu_k_vec);
totalDuration = numPeriods * Ts;

%----------------------------------------------------------------------
% Event Generation Across Multiple Periods
%----------------------------------------------------------------------
% Accumulators for absolute event times across all periods
detectedAll = [];
missedAll = [];
deadtimeBlocks = [];
lastDetectionOffset = -tau - 1; % Negative offset => detector recovered before t=0

for i = 1:numPeriods
    currentMu_k = mu_k_vec(i);
    
    % Simulate events within this period (times are relative to period start)
    [detectedTime, missedTime, ~] = SPAD_detection_complex(currentMu_k, Ts, tau, lastDetectionOffset);
    startTime = (i - 1) * Ts;
    
    % Convert per-period times to absolute times for accumulation
    detectedAbs = detectedTime + startTime;
    missedAbs = missedTime + startTime;
    
    detectedAll = [detectedAll, detectedAbs];
    missedAll = [missedAll, missedAbs];
    
    % Each detection creates a dead-time interval clipped to total duration
    deadtimeBlocks = [deadtimeBlocks; [detectedAbs', min(detectedAbs' + tau, totalDuration)]];
    
    % Update dead time carryover for next period
    if isempty(detectedTime)
        % No detections: shift offset by one period, but keep it capped
        lastDetectionOffset = max(lastDetectionOffset - Ts, -tau - 1);
    else
        % Carry the last detection time into the next period
        latestDetectedTime = startTime + max(detectedTime);
        lastDetectionOffset = latestDetectedTime - (i * Ts);
    end
end

%----------------------------------------------------------------------
% Create Publication-Quality Figure
%----------------------------------------------------------------------
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 7.4 3.2]);
ax = axes('Parent', fig);
hold(ax, 'on');

% Color scheme
periodColor = [0.88 0.92 0.97];      % Light blue for periods
recoveryColor = [0.55 0.70 0.82];    % Darker blue for dead time
detColor = [0.05 0.45 0.75];         % Dark blue for detections
missColor = [0.85 0.20 0.25];        % Red for missed events
timelineColor = [0.35 0.35 0.35];    % Gray for timeline
yMin = -0.17;
yMax = 0.17;

% Draw period backgrounds
for k = 0:numPeriods-1
    patch(ax, [k*Ts, (k+1)*Ts, (k+1)*Ts, k*Ts], [yMin, yMin, yMax, yMax], ...
        periodColor, 'EdgeColor', 'none', 'FaceAlpha', 0.45, 'HandleVisibility', 'off');
end

% Draw dead time regions
for idx = 1:size(deadtimeBlocks, 1)
    patch(ax, [deadtimeBlocks(idx, 1), deadtimeBlocks(idx, 2), deadtimeBlocks(idx, 2), deadtimeBlocks(idx, 1)], ...
        [yMin, yMin, yMax, yMax], recoveryColor, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

% Draw timeline
plot(ax, [0, totalDuration], [0, 0], '--', 'Color', timelineColor, 'LineWidth', 1.1, 'HandleVisibility', 'off');

% Plot missed events (crosses)
scatter(ax, missedAll, zeros(size(missedAll)), 700, missColor, 'LineWidth', 3, ...
    'Marker', 'x', 'HandleVisibility', 'off');

% Plot detected events (filled circles)
scatter(ax, detectedAll, zeros(size(detectedAll)), 700, detColor, 'filled', ...
    'MarkerFaceAlpha', 1, 'MarkerEdgeColor', 'w', 'LineWidth', 2, 'HandleVisibility', 'off');

% Draw period boundaries
for k = 0:numPeriods
    xline(ax, k * Ts, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8, 'HandleVisibility', 'off');
end

%----------------------------------------------------------------------
% Figure Formatting
%----------------------------------------------------------------------
xlabel(ax, 'Time (ns)', 'FontWeight', 'bold');

set(ax, 'YTick', [], 'YLim', [yMin, yMax], 'XLim', [0, totalDuration], ...
    'FontName', 'Arial', 'FontSize', 11, 'LineWidth', 1, 'Box', 'on', 'Layer', 'top');

grid(ax, 'off');

%----------------------------------------------------------------------
% Export to PDF (Vector Format)
%----------------------------------------------------------------------
drawnow;
exportgraphics(fig, 'SPAD_detection_publication.pdf', 'ContentType', 'vector');

fprintf('=== Simulation Complete ===\n');
fprintf('Total periods: %d\n', numPeriods);
fprintf('Detected events: %d\n', length(detectedAll));
fprintf('Missed events: %d\n', length(missedAll));
fprintf('Output saved to: SPAD_detection_publication.pdf\n');
