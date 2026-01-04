%==========================================================================
% SiPM 30020 Full-Array Monte Carlo Simulation
%==========================================================================
% Description:
%   Simulates SiPM (Silicon Photomultiplier) response using full-array
%   Monte Carlo where each SPAD (Single Photon Avalanche Diode) is
%   explicitly simulated. Models photon detection efficiency, dead time,
%   and saturation effects. Compares simulation results with analytic
%   model. Note: Very slow for large arrays.
%
% Device Model:
%   SiPM 30020 (14410 SPADs, 3.07×3.07 mm² active area)
%
% Performance Metrics:
%   - Count rate vs. irradiance
%   - Bias current vs. irradiance
%   - Comparison with analytic saturation model
%
% Author: Akatsuki Sky
% Date: 2026-1-4
%==========================================================================

clc; clear;

%----------------------------------------------------------------------
% Device Parameters (SiPM 30020)
%----------------------------------------------------------------------
opticalIrradiance = logspace(-4, 1, 20); % incident irradiance at the SiPM surface (W/m^2)
PDE               = 0.38;                % photon detection efficiency (probability a photon triggers a cell)
lambda            = 405e-9;             % illumination wavelength used for photon energy (m)
c                 = 3e8;                % speed of light for photon energy (m/s)
h                 = 6.626e-34;          % Planck constant for photon energy (J*s)
Ts                = 15e-9;              % counting window duration per symbol (s)
deadTime          = 60e-9;              % SPAD recovery/dead time after a detection (s)
SiPMArea          = (3.07e-3)^2;        % total active area of the SiPM (m^2)
Ncells            = 14410;              % number of microcells (SPADs) in the array
Q_cell            = 0.17e-12;           % charge released per fired microcell (C)

symbols = 100; % Number of counting windows per trial (kept small for runtime)
trials  = 10;  % Number of Monte Carlo trials per irradiance point (kept small)

%----------------------------------------------------------------------
% Analytic Model (Dead Time Saturation)
%----------------------------------------------------------------------
% Photon energy and responsivity
Ep      = h * c / lambda;
A_pixel = SiPMArea / Ncells;
alpha   = (PDE * A_pixel) / Ep;

% Analytic bias current (Eq. 5 in paper)
I_analytic = (Q_cell * Ncells * alpha .* opticalIrradiance) ./ (1 + alpha * deadTime .* opticalIrradiance);

% Analytic count rate (Eq. 3 in paper)
opticalPower = opticalIrradiance .* SiPMArea;
photonEnergy = c * h / lambda;
R_in         = opticalPower ./ photonEnergy;
R_ideal      = R_in .* PDE;
R_analytic   = R_ideal ./ (1 + (R_ideal .* deadTime) ./ Ncells);
R_analytic_Gcps = R_analytic / 1e9;

%----------------------------------------------------------------------
% Full-Array Monte Carlo Simulation
%----------------------------------------------------------------------
aveArrive        = (opticalPower .* Ts) * PDE / photonEnergy;
aveArrivePerSPAD = aveArrive ./ Ncells; % expected arrivals per SPAD per window
aveCount = zeros(size(aveArrivePerSPAD));

fprintf('Full-array MC (slow). %d points, %d trials, %d periods.\n', ...
    numel(opticalIrradiance), trials, symbols);

for q = 1:numel(aveArrivePerSPAD)
    start_rate = aveArrivePerSPAD(q);
    trialMean  = zeros(trials, 1);
    
    for j = 1:trials
        startTime = zeros(Ncells, 1); % Per-SPAD carry-over dead-time offset
        counts = zeros(symbols, 1);
        
        for i = 1:symbols
            totalCount = 0;
            % Simulate each SPAD individually
            for spadIdx = 1:Ncells
                [countSpad, startTime(spadIdx)] = Count_func(start_rate, Ts, deadTime, startTime(spadIdx));
                totalCount = totalCount + countSpad;
            end
            counts(i) = totalCount;
        end
        trialMean(j) = mean(counts);
    end
    
    aveCount(q) = mean(trialMean);
    fprintf('Point %d/%d done.\n', q, numel(aveArrivePerSPAD));
end

% Convert average counts to rate and current
R_sim_Gcps = aveCount / (1e9 * Ts); % Count rate in Gcps
I_mc = (Q_cell * aveCount) / Ts;    % Bias current in Amperes

%----------------------------------------------------------------------
% Plot 1: Count Rate vs. Irradiance
%----------------------------------------------------------------------
figure('Name','SiPM 30020 Count Rate (Full)','Color','w','Position',[100 100 880 560]);
loglog(opticalIrradiance, R_analytic_Gcps, 'Color',[0.80 0.10 0.10], ...
    'LineWidth',2.2, 'LineStyle','-', 'DisplayName','Analytic (Eq. 3)');
hold on; grid on;
loglog(opticalIrradiance, R_sim_Gcps, 'Color',[0.00 0.55 0.20], ...
    'LineWidth',2.0, 'LineStyle','--', 'Marker','o', 'MarkerSize',7, ...
    'MarkerFaceColor',[0.85 1.00 0.90], 'DisplayName','Full-array MC');
loglog(opticalIrradiance, R_ideal/1e9, 'Color',[0.30 0.30 0.30], ...
    'LineWidth',1.6, 'LineStyle',':', 'DisplayName','Ideal (Linear)');

C_max_Gcps = (Ncells / deadTime)/1e9;
yline(C_max_Gcps, ':', 'Color',[0.5 0.5 0.5], 'HandleVisibility','off');

xlabel('Irradiance (W/m^2)', 'FontName','Times','FontSize',14);
ylabel('Count Rate (Gcps)', 'FontName','Times','FontSize',14);
title('SiPM 30020 Count Rate: Full MC vs Analytic', 'FontName','Times','FontSize',16);
legend('Location', 'best', 'FontSize', 11);
set(gca,'FontName','Times','FontSize',14);

%----------------------------------------------------------------------
% Plot 2: Bias Current vs. Irradiance
%----------------------------------------------------------------------
figure('Name','SiPM 30020 Bias Current (Full)','Color','w','Position',[120 120 880 560]);
loglog(opticalIrradiance, I_analytic, 'Color',[0.80 0.10 0.10], ...
    'LineWidth',2.2, 'LineStyle','-', 'DisplayName','Analytic (Eq. 5)');
hold on; grid on;
loglog(opticalIrradiance, I_mc, 'Color',[0.00 0.55 0.20], ...
    'LineWidth',2.0, 'LineStyle','--', 'Marker','o', 'MarkerSize',7, ...
    'MarkerFaceColor',[0.85 1.00 0.90], 'DisplayName','Full-array MC');

xlabel('Irradiance (W/m^2)', 'FontName','Times','FontSize',14);
ylabel('Bias Current (A)', 'FontName','Times','FontSize',14);
title('SiPM 30020 Bias Current: Full MC vs Analytic', 'FontName','Times','FontSize',16);
legend('Location', 'best', 'FontSize', 11);
set(gca,'FontName','Times','FontSize',14);
