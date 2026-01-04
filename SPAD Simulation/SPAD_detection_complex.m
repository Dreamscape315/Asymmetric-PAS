%==========================================================================
% SPAD Detection Simulation with Dead Time
%==========================================================================
% Description:
%   Simulates Single Photon Avalanche Diode (SPAD) detection events
%   within one counting period. Models non-paralyzable dead time where
%   the detector is unresponsive after each detection. Uses Poisson
%   process for photon arrival modeling.
%
% Inputs:
%   mu_k - Expected photon count in one period (Poisson parameter)
%   Ts - Period length (ns)
%   tau - Detector dead time after each detection (ns)
%   lastDetectedOffset - Last detection time relative to period start
%                        (negative if detection was in previous period)
%
% Outputs:
%   detectedTimes - Array of detected event times within [0, Ts)
%   missedTimes - Array of event times missed due to dead time
%   numDetected - Total number of detected events
%
% Examples:
%   [det, miss, num] = SPAD_detection_complex(10, 180, 45, -50);
%
% Author: Akatsuki Sky
% Date: 2026-1-1
%==========================================================================

function [detectedTimes, missedTimes, numDetected] = SPAD_detection_complex(mu_k, Ts, tau, lastDetectedOffset)
    %----------------------------------------------------------------------
    % Initialization
    %----------------------------------------------------------------------
    detectedTimes = [];
    missedTimes = [];
    currentTime = 0;
    lambda = mu_k / Ts; % Poisson process rate (events per ns)
    
    %----------------------------------------------------------------------
    % Event Generation Loop (Poisson Process)
    %----------------------------------------------------------------------
    while true
        % Generate next photon arrival time using exponential distribution
        u = rand();
        nextArrivalTime = currentTime - (1 / lambda) * log(1 - u);
        
        % Check if event occurs within current period
        if nextArrivalTime >= Ts
            break;
        end
        
        currentTime = nextArrivalTime;
        timeSinceLastDetection = currentTime - lastDetectedOffset;
        
        % Check detector state
        if timeSinceLastDetection <= tau
            % Event arrives during dead time: missed detection
            missedTimes = [missedTimes, currentTime];
        else
            % Detector has recovered: record detection and enter new dead time
            detectedTimes = [detectedTimes, currentTime];
            lastDetectedOffset = currentTime;
        end
    end
    
    numDetected = length(detectedTimes);
end
