%% project1_link_budget.m
%  Project 1: Starlink LEO Link Budget + Patch Array Antenna Simple Design
%  Course  : Defense Antenna Engineering
%  Student : 202550592 오선지
%  Theory  : Lecture Note 7.1 Ch.5, 9.1 / Balanis Antenna Theory 3rd Ed.
%
%  Run from the script directory.  Outputs:
%    project1_linkbudget.png   — link budget bar chart
%    project1_array.png        — patch array layout + gain curve

clear; clc; close all;

%% ═══ Physical Constants ══════════════════════════════════════════════════
c  = 299792458;      % speed of light [m/s]
kB = 1.3807e-23;     % Boltzmann constant [J/K]
T0 = 290.0;          % reference temperature [K]

%% ═══ MCS Table (802.11ac, 40 MHz BW, 1SS, GI=800ns) ════════════════════
%  From the provided Excel modem spec sheet  (1-based index)
MCS_mod  = {'BPSK','QPSK','QPSK','16-QAM','16-QAM', ...
            '64-QAM','64-QAM','64-QAM','256-QAM','256-QAM'};
MCS_cod  = {'1/2','1/2','3/4','1/2','3/4','2/3','3/4','5/6','3/4','5/6'};
MCS_rate = [13.5,27.0,40.5,54.0,81.0,108.0,121.5,135.0,162.0,180.0]; % Mbps
MCS_snr  = [5,8,12,14,18,21,23,28,32,34];                             % dB

%% ═══ Starlink Gen1 System Parameters ════════════════════════════════════
%  Source: FCC STA filings, ITU coordination documents

% Orbital
ALT_km   = 550.0;
ELEV_deg = 40.0;
SLANT_km = ALT_km / sind(ELEV_deg);   % slant range [km]
SLANT_m  = SLANT_km * 1e3;

% ── Downlink: Satellite → User Terminal (Ku-band 12 GHz) ──────────────
DL_f_GHz  = 12.0;
DL_f      = DL_f_GHz * 1e9;
DL_lam    = c / DL_f;
DL_BW_MHz = 40.0;
DL_BW     = DL_BW_MHz * 1e6;
DL_MCS    = 8;   % 1-based (MCS7: 64-QAM 5/6, req SNR = 28 dB)

SAT_EIRP_dBW   = 57.0;                 % [dBW] from FCC STA
SAT_TX_P_dBm   = 43.0;                 % ~20 W
SAT_TX_G_dBi   = SAT_EIRP_dBW - (SAT_TX_P_dBm - 30);

UT_DL_RX_G_dBi = 38.0;                 % UT patch array target
UT_DL_NF_dB    = 1.5;                  % LNA noise figure

% ── Uplink: User Terminal → Satellite (Ku-band 14.25 GHz) ─────────────
UL_f_GHz  = 14.25;
UL_f      = UL_f_GHz * 1e9;
UL_lam    = c / UL_f;
UL_BW_MHz = 40.0;
UL_BW     = UL_BW_MHz * 1e6;
UL_MCS    = 4;   % 1-based (MCS3: 16-QAM 1/2, req SNR = 14 dB)

UT_UL_TX_P_dBm = 33.0;                 % ~2 W
UT_UL_TX_G_dBi = 36.0;

SAT_GT_dBK   =  3.0;                   % [dB/K]
T_SAT        = 700.0;                  % [K] (facing warm Earth)
SAT_RX_G_dBi = SAT_GT_dBK + 10*log10(T_SAT);
SAT_NF_dB    = 10*log10(1 + (T_SAT - T0)/T0);

% ── Additional Losses ─────────────────────────────────────────────────
LOSS_BEAM_dB  = 1.0;   % beam-pointing mismatch
LOSS_POL_dB   = 0.5;   % polarization mismatch (PLF)
LOSS_IMP_dB   = 0.2;   % impedance mismatch
LOSS_CABLE_dB = 0.3;   % cable / connector
LOSS_FADE_dB  = 3.0;   % fading / atmospheric margin
TOTAL_LOSS_dB = LOSS_BEAM_dB + LOSS_POL_dB + ...
                LOSS_IMP_dB  + LOSS_CABLE_dB + LOSS_FADE_dB;

%% ═══ Link Budget ══════════════════════════════════════════════════════════
%  Friis (dB form):   P_rx = P_tx + G_tx − FSPL + G_rx
%  FSPL [dB]      = 20·log10(4πRf/c)
%  N_eff [dBm]    = 10·log10(k·T0·B) + 30 + NF
%  SNR_final [dB] = P_rx − N_eff − ΣLoss

DL_fspl   = 20*log10(4*pi*SLANT_m*DL_f/c);
UL_fspl   = 20*log10(4*pi*SLANT_m*UL_f/c);

DL_rx_pwr = SAT_TX_P_dBm + SAT_TX_G_dBi - DL_fspl + UT_DL_RX_G_dBi;
UL_rx_pwr = UT_UL_TX_P_dBm + UT_UL_TX_G_dBi - UL_fspl + SAT_RX_G_dBi;

DL_noise  = 10*log10(kB*T0*DL_BW) + 30 + UT_DL_NF_dB;
UL_noise  = 10*log10(kB*T0*UL_BW) + 30 + SAT_NF_dB;

DL_snr_raw = DL_rx_pwr - DL_noise;
DL_snr_fin = DL_snr_raw - TOTAL_LOSS_dB;
DL_margin  = DL_snr_fin - MCS_snr(DL_MCS);

UL_snr_raw = UL_rx_pwr - UL_noise;
UL_snr_fin = UL_snr_raw - TOTAL_LOSS_dB;
UL_margin  = UL_snr_fin - MCS_snr(UL_MCS);

%% ── Print ─────────────────────────────────────────────────────────────
fprintf('\n%s\n  DOWNLINK  (Satellite → UT, %.0f GHz)\n%s\n', ...
        repmat('=',1,60), DL_f_GHz, repmat('=',1,60));
fprintf('  FSPL              = %7.2f dB\n', DL_fspl);
fprintf('  Received Power    = %7.2f dBm\n', DL_rx_pwr);
fprintf('  Noise Floor       = %7.2f dBm\n', DL_noise);
fprintf('  SNR (raw)         = %7.2f dB\n', DL_snr_raw);
fprintf('  Total Loss        = %7.2f dB\n', TOTAL_LOSS_dB);
fprintf('  Final SNR         = %7.2f dB\n', DL_snr_fin);
fprintf('  Required SNR      = %7.0f dB   [MCS%d %s %s, %.1f Mbps]\n', ...
        MCS_snr(DL_MCS), DL_MCS-1, MCS_mod{DL_MCS}, MCS_cod{DL_MCS}, MCS_rate(DL_MCS));
if DL_margin >= 0
    fprintf('  Margin            = %7.2f dB   --> PASS\n', DL_margin);
else
    fprintf('  Margin            = %7.2f dB   --> FAIL\n', DL_margin);
end

fprintf('\n%s\n  UPLINK  (UT → Satellite, %.2f GHz)\n%s\n', ...
        repmat('=',1,60), UL_f_GHz, repmat('=',1,60));
fprintf('  FSPL              = %7.2f dB\n', UL_fspl);
fprintf('  Received Power    = %7.2f dBm\n', UL_rx_pwr);
fprintf('  Noise Floor       = %7.2f dBm\n', UL_noise);
fprintf('  SNR (raw)         = %7.2f dB\n', UL_snr_raw);
fprintf('  Total Loss        = %7.2f dB\n', TOTAL_LOSS_dB);
fprintf('  Final SNR         = %7.2f dB\n', UL_snr_fin);
fprintf('  Required SNR      = %7.0f dB   [MCS%d %s %s, %.1f Mbps]\n', ...
        MCS_snr(UL_MCS), UL_MCS-1, MCS_mod{UL_MCS}, MCS_cod{UL_MCS}, MCS_rate(UL_MCS));
if UL_margin >= 0
    fprintf('  Margin            = %7.2f dB   --> PASS\n', UL_margin);
else
    fprintf('  Margin            = %7.2f dB   --> FAIL\n', UL_margin);
end

%% ═══ Patch Array Antenna Simple Design ════════════════════════════════════
%  Substrate: Rogers RT/duroid 5880  (εr=2.2, h=1.57 mm)
%  Balanis "Antenna Theory" Ch.14:
%   W      = (c/2f)·√(2/(εr+1))
%   εr_eff = (εr+1)/2 + (εr-1)/2·(1+12h/W)^(-½)
%   ΔL     = 0.412h·[(εr_eff+0.3)(W/h+0.264)] / [(εr_eff-0.258)(W/h+0.8)]
%   L      = c/(2f√εr_eff) − 2ΔL
%   G_array= G_el + 10·log10(Nx·Ny)

er   = 2.2;
h    = 1.57e-3;
G_el = 7.0;

ant_DL = patch_array_design(DL_f, UT_DL_RX_G_dBi, er, h, G_el, c, 'DL Rx @ 12.0 GHz');
ant_UL = patch_array_design(UL_f, UT_UL_TX_G_dBi, er, h, G_el, c, 'UL Tx @ 14.25 GHz');

%% ═══ Figure 1: Link Budget Summary ════════════════════════════════════════
fig1 = figure('Name','Link Budget Summary','Position',[50 50 900 500]);

labels = {'DL SNR (raw)','DL Total Loss','DL Final SNR','DL Req SNR','DL Margin', ...
          'UL SNR (raw)','UL Total Loss','UL Final SNR','UL Req SNR','UL Margin'};
vals   = [DL_snr_raw, TOTAL_LOSS_dB, DL_snr_fin, MCS_snr(DL_MCS), DL_margin, ...
          UL_snr_raw, TOTAL_LOSS_dB, UL_snr_fin, MCS_snr(UL_MCS), UL_margin];

% Color each bar
clrs = [0.27 0.51 0.71;  % DL SNR
        0.93 0.40 0.40;  % DL Loss
        0.27 0.51 0.71;  % DL Final
        1.00 0.65 0.00;  % DL Req
        0.00 0.65 0.30;  % DL Margin
        0.39 0.60 0.93;  % UL SNR
        0.93 0.40 0.40;  % UL Loss
        0.39 0.60 0.93;  % UL Final
        1.00 0.65 0.00;  % UL Req
        0.00 0.65 0.30]; % UL Margin
if DL_margin < 0, clrs(5,:) = [0.85 0.15 0.15]; end
if UL_margin < 0, clrs(10,:) = [0.85 0.15 0.15]; end

b = bar(vals, 'FaceColor','flat');
b.CData = clrs;
set(gca, 'XTick',1:10, 'XTickLabel',labels, 'XTickLabelRotation',35, 'FontSize',9);
ylabel('Value (dB)', 'FontSize',11);
title(sprintf('Starlink Link Budget  —  DL %.0f GHz | UL %.2f GHz', ...
              DL_f_GHz, UL_f_GHz), 'FontSize',13);
hold on;
plot([0.5 10.5],[0 0],'k-','LineWidth',1.0);
hold off;
grid on;
% Value labels above each bar
for i = 1:numel(vals)
    if vals(i) >= 0
        text(i, vals(i)+0.6, sprintf('%.1f',vals(i)), ...
             'HorizontalAlignment','center','FontSize',8);
    else
        text(i, vals(i)-1.5, sprintf('%.1f',vals(i)), ...
             'HorizontalAlignment','center','FontSize',8);
    end
end

print(fig1, '-dpng', '-r150', 'project1_linkbudget.png');
fprintf('\n[Saved]  project1_linkbudget.png\n');

%% ═══ Figure 2: Patch Array Layout + Gain Curve ═════════════════════════════
fig2 = figure('Name','Patch Array Design','Position',[100 100 1100 480]);

% ── (a) Array layout (DL 12 GHz, show up to 10×10 sub-array) ─────────────
subplot(1,2,1);
show = min([10, ant_DL.Nx, ant_DL.Ny]);
hold on;
for ix = 0:show-1
    for iy = 0:show-1
        x0 = ix * ant_DL.d_mm;
        y0 = iy * ant_DL.d_mm;
        rectangle('Position',[x0 - ant_DL.W_mm/2,  y0 - ant_DL.L_mm/2, ...
                               ant_DL.W_mm,          ant_DL.L_mm], ...
                  'EdgeColor',[0 0 0.5], 'FaceColor',[0.53 0.81 0.98], ...
                  'LineWidth',0.6);
    end
end
hold off;
axis equal; grid on;
xlabel('x  (mm)','FontSize',11);
ylabel('y  (mm)','FontSize',11);
title(sprintf(['DL Rx Patch Array  (12 GHz)\n' ...
               'Full: %d×%d = %d elements  —  %.0f×%.0f cm\n' ...
               'Showing %d×%d sub-array'], ...
    ant_DL.Nx, ant_DL.Ny, ant_DL.Nx*ant_DL.Ny, ...
    ant_DL.Ax_cm, ant_DL.Ay_cm, show, show), 'FontSize',10);

% Patch dimension annotation
text(show*ant_DL.d_mm*0.5, -ant_DL.d_mm*0.7, ...
     sprintf('W=%.2f mm,  L=%.2f mm,  d=%.2f mm', ...
             ant_DL.W_mm, ant_DL.L_mm, ant_DL.d_mm), ...
     'HorizontalAlignment','center','FontSize',8,'Color',[0 0 0.5]);

% ── (b) Gain vs Number of Elements ────────────────────────────────────────
subplot(1,2,2);
N_vec = 1:5000;
G_vec = G_el + 10*log10(N_vec);
semilogx(N_vec, G_vec, 'b-', 'LineWidth',1.8); hold on;

% Target gain lines
plot([1 5000],[UT_DL_RX_G_dBi UT_DL_RX_G_dBi],'g--','LineWidth',1.5);
plot([1 5000],[UT_UL_TX_G_dBi UT_UL_TX_G_dBi],'r--','LineWidth',1.5);

% Mark design points
N_DL = ant_DL.Nx * ant_DL.Ny;
N_UL = ant_UL.Nx * ant_UL.Ny;
plot(N_DL, ant_DL.G_dBi, 'go', 'MarkerSize',10, 'LineWidth',2);
plot(N_UL, ant_UL.G_dBi, 'rs', 'MarkerSize',10, 'LineWidth',2);
hold off;

xlabel('Number of Elements  (N_x \times N_y)','FontSize',11);
ylabel('Array Gain  (dBi)','FontSize',11);
title('Patch Array Gain vs Number of Elements','FontSize',12);
legend('G_{el}=7 dBi', ...
       sprintf('DL target %.0f dBi',UT_DL_RX_G_dBi), ...
       sprintf('UL target %.0f dBi',UT_UL_TX_G_dBi), ...
       sprintf('DL: %d el → %.1f dBi',N_DL,ant_DL.G_dBi), ...
       sprintf('UL: %d el → %.1f dBi',N_UL,ant_UL.G_dBi), ...
       'Location','northwest','FontSize',9);
xlim([1 5000]); grid on;

print(fig2, '-dpng', '-r150', 'project1_array.png');
fprintf('[Saved]  project1_array.png\n');

fprintf('\n%s\n  FINAL SUMMARY\n%s\n', repmat('=',1,60), repmat('=',1,60));
fprintf('  DL: FSPL=%5.1fdB  SNR=%5.1fdB  Margin=%+5.1fdB  %s\n', ...
        DL_fspl, DL_snr_fin, DL_margin, sel('PASS','FAIL',DL_margin>=0));
fprintf('  UL: FSPL=%5.1fdB  SNR=%5.1fdB  Margin=%+5.1fdB  %s\n', ...
        UL_fspl, UL_snr_fin, UL_margin, sel('PASS','FAIL',UL_margin>=0));
fprintf('  DL Array: %dx%d=%d el,  W=%.2fmm  L=%.2fmm,  G=%.1fdBi\n', ...
        ant_DL.Nx,ant_DL.Ny,N_DL,ant_DL.W_mm,ant_DL.L_mm,ant_DL.G_dBi);
fprintf('  UL Array: %dx%d=%d el,  W=%.2fmm  L=%.2fmm,  G=%.1fdBi\n', ...
        ant_UL.Nx,ant_UL.Ny,N_UL,ant_UL.W_mm,ant_UL.L_mm,ant_UL.G_dBi);
fprintf('%s\n', repmat('=',1,60));


%% ═══ Local Functions ══════════════════════════════════════════════════════

function ant = patch_array_design(freq, G_tgt_dBi, er, h, G_el_dBi, c, label)
% Simple patch array design  (Balanis Ch.14 / Lecture Note 7.1 Ch.5)
    lam = c / freq;

    % ① Patch width  W  [Eq.14-7]
    W = (c / (2*freq)) * sqrt(2 / (er + 1));

    % ② Effective permittivity  εr_eff  [Eq.14-1]
    er_eff = (er+1)/2 + (er-1)/2 * (1 + 12*h/W)^(-0.5);

    % ③ End extension  ΔL  [Eq.14-2]
    dL = 0.412 * h * ((er_eff+0.3)*(W/h+0.264)) / ...
                     ((er_eff-0.258)*(W/h+0.8));

    % ④ Physical length  L
    L = lam / (2*sqrt(er_eff)) - 2*dL;

    % ⑤ Element spacing  d = λ/2
    d = lam / 2;

    % ⑥ Number of elements  (square array, rounded up to even)
    N_side = ceil(sqrt(10^((G_tgt_dBi - G_el_dBi)/10)));
    if mod(N_side,2), N_side = N_side + 1; end
    Nx = N_side;  Ny = N_side;

    % ⑦ Array gain
    G_array_dBi = G_el_dBi + 10*log10(Nx * Ny);

    fprintf('\n--- %s ---\n', label);
    fprintf('  W = %.3f mm,  εr_eff = %.4f,  ΔL = %.3f mm,  L = %.3f mm\n', ...
            W*1e3, er_eff, dL*1e3, L*1e3);
    fprintf('  d = %.2f mm (λ/2),  Array: %dx%d = %d el\n', d*1e3, Nx, Ny, Nx*Ny);
    fprintf('  Aperture: %.1f × %.1f cm,  Gain = %.2f dBi\n', ...
            Nx*d*100, Ny*d*100, G_array_dBi);

    ant.W_mm   = W*1e3;
    ant.L_mm   = L*1e3;
    ant.er_eff = er_eff;
    ant.dL_mm  = dL*1e3;
    ant.d_mm   = d*1e3;
    ant.Nx     = Nx;
    ant.Ny     = Ny;
    ant.G_dBi  = G_array_dBi;
    ant.Ax_cm  = Nx*d*100;
    ant.Ay_cm  = Ny*d*100;
end

function s = sel(a, b, cond)
    if cond, s = a; else, s = b; end
end
