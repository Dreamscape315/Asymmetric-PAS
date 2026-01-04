%==========================================================================
% SPAD Photon Counting with Non-Paralyzable Dead Time
%==========================================================================
% Description:
%   Simulates photon detection by a SPAD (Single Photon Avalanche Diode)
%   within one counting period. Models non-paralyzable dead time where
%   the detector is unresponsive after each detection for a fixed duration.
%   Uses Poisson process for photon arrivals.
%
% Inputs:
%   average_photon_number - Expected photon count in one period (lambda)
%   counting_period - Duration of counting window
%   tau - Dead time after each detection
%   start_time - Time offset from previous period's dead time
%
% Outputs:
%   numb_count - Number of detected photons in this period
%   Interference_time - Residual dead time carried to next period
%
% Examples:
%   [count, carryover] = Count_func(10, 100, 5, 0);
%
% Author: Cuiwei HE (Original)
% Editor: Akatsuki Sky
% Date: 2026-1-4
%==========================================================================

function [numb_count,Interference_time]= Count_func(average_photon_number,counting_period,tau,start_time)

    %----------------------------------------------------------------------
    % Early Exit Conditions
    %----------------------------------------------------------------------
    % Case 1: Still in dead time from previous period
    if start_time>=counting_period
        numb_count=0;
        Interference_time=start_time-counting_period;
        return;
    end
    
    % Case 2: No photons expected
    if average_photon_number<=0
        numb_count=0;
        Interference_time=0;
        return;
    end
    
    %----------------------------------------------------------------------
    % Photon Detection Simulation (Poisson Process)
    %----------------------------------------------------------------------
    rate_parameter=average_photon_number/counting_period;
    time=start_time;
    last_detected=-Inf;
    numb_count=0;
    
    % Simulate photon arrivals using exponential inter-arrival times
    while true
        % Generate next photon arrival time
        time=time+(-log(1-rand)/rate_parameter);
        
        % Check if photon arrives within counting period
        if time>=counting_period
            break;
        end
        
        % Check if detector has recovered from dead time
        if time-last_detected>=tau
            numb_count=numb_count+1;
            last_detected=time;
        end
    end
    
    %----------------------------------------------------------------------
    % Calculate Dead Time Carryover to Next Period
    %----------------------------------------------------------------------
    if numb_count>0
        T_gap=tau-(counting_period-last_detected);
        if T_gap>0
            Interference_time=T_gap;
        else
            Interference_time=0;
        end
    else
        Interference_time=0;
    end
end
