%% project2_far_field.m
%  Project 2: Far-Field Pattern from Planar Near-Field Measurements
%  Course  : Defense Antenna Engineering
%  Student : 202550592 오선지
%  AUT     : Standard Gain Horn (SGH) @ 5.3 GHz
%
%  Theory  : Lecture Note 11.1 / Balanis "Antenna Theory" Ch.12.9
%  Method  :
%    Step 1 — Back-Projection  : aperture field recovery
%    Step 2 — Plane Wave Spectrum (PWS) via 2-D IFFT
%    Step 3 — NF→FF on Ewald sphere
%             E_θ ∝ cosθ[cosφ·fx + sinφ·fy]
%             E_φ ∝      [−sinφ·fx + cosφ·fy]
%    Step 4 — 1-D cuts: E-plane (φ=0°), H-plane (φ=90°)
%
%  Outputs (saved in working directory):
%    project2_backproj.png     — back-projection surface
%    project2_farfield.png     — 1-D far-field pattern (Cartesian + Polar)
%    project2_pws.png          — plane wave spectrum
%
%  Run from:  .../pr2_files/   (where SGH_FullPlanarScan.mat lives)
%  OR set data_dir below.

clear; clc; close all;

%% ─── Path & Data ────────────────────────────────────────────────────────
script_dir = fileparts(mfilename('fullpath'));
data_dir   = fullfile(script_dir, 'pr2_files');
addpath(data_dir);

load(fullfile(data_dir, 'SGH_FullPlanarScan.mat'));

%% ─── Measurement Parameters ─────────────────────────────────────────────
c   = 299792458;
f   = 5.3e9;        % frequency [Hz]
lam = c / f;        % wavelength [m]  (56.56 mm)
k0  = 2*pi / lam;   % wavenumber [rad/m]
z_0 = 0.088;        % probe–AUT distance [m]
dx  = 0.022;        % x scan step [m]
dy  = 0.022;        % y scan step [m]

E_meas_x = comp_x_measurement(:,:,3);   % (68×36) complex, @ 5.3 GHz
E_meas_y = comp_y_measurement(:,:,3);
[M, N]   = size(E_meas_x);
MI = 10*M;   NI = 10*N;                 % zero-padded FFT size

fprintf('f = %.2f GHz | lam = %.2f mm | k0 = %.2f rad/m\n', f/1e9, lam*1e3, k0);
fprintf('Scan grid   : %d x %d  (dx=%.0f mm)\n', M, N, dx*1e3);
fprintf('FFT size    : %d x %d  (10x zero-padded)\n', MI, NI);

%% ═══════════════════════════════════════════════════════════════════════════
%%  STEP 1 — Back-Projection  (Balanis Eq. 12-86)
%%  Uses provided BackProjection_PlanarNearField.m
%% ═══════════════════════════════════════════════════════════════════════════
[Ex_bp, Ey_bp, Ez_bp, x_ax, y_ax] = ...
    BackProjection_PlanarNearField(E_meas_x, E_meas_y, z_0, f, dx, dy);

Ex_bp = abs(fliplr(Ex_bp'));   % (N×M) after transpose + flip → view from front
Ey_bp = abs(fliplr(Ey_bp'));
Ez_bp = abs(fliplr(Ez_bp'));
U_bp  = Ex_bp.^2 + Ey_bp.^2 + Ez_bp.^2;
U_dB  = 10*log10(U_bp ./ max(U_bp(:)) + 1e-10);

%% ═══════════════════════════════════════════════════════════════════════════
%%  STEP 2 — Plane Wave Spectrum  (Balanis Eq. 12-85)
%%  e^(+jωt) convention → use ifft2 (not fft2)
%%    fx(kx,ky) = IFFT2D { Ex_meas(x,y,z0) }
%% ═══════════════════════════════════════════════════════════════════════════
fx = ifftshift(ifft2(E_meas_x, MI, NI));   % (MI×NI)
fy = ifftshift(ifft2(E_meas_y, MI, NI));

m  = (-MI/2 : MI/2-1);          % row  frequency indices
n  = (-NI/2 : NI/2-1);          % col  frequency indices
kx = 2*pi*m / (MI*dx);          % [rad/m], 1×MI
ky = 2*pi*n / (NI*dy);          % [rad/m], 1×NI

fprintf('PWS kx range: +/-%.1f rad/m  (k0=%.2f rad/m)\n', max(abs(kx)), k0);

%% ═══════════════════════════════════════════════════════════════════════════
%%  STEP 3 — NF→FF Transformation  (Lecture Note 11.1 / Balanis §12.9)
%%
%%  Far-field via plane wave spectrum (sampled on Ewald sphere):
%%    E_θ(r,θ,φ) ∝ cosθ · [cosφ·fx(kx,ky) + sinφ·fy(kx,ky)]
%%    E_φ(r,θ,φ) ∝        [−sinφ·fx(kx,ky) + cosφ·fy(kx,ky)]
%%    where  kx = k0·sinθ·cosφ,  ky = k0·sinθ·sinφ
%%
%%  E-plane (φ=0°) : kx = k0·sinθ,  ky = 0
%%  H-plane (φ=90°): kx = 0,  ky = k0·sinθ
%% ═══════════════════════════════════════════════════════════════════════════
theta_deg = linspace(-90, 90, 1801);
theta_rad = deg2rad(theta_deg);
sin_th    = sin(theta_rad);
cos_th    = cos(theta_rad);

% DC indices (1-based)
ky0 = NI/2 + 1;
kx0 = MI/2 + 1;

% ── E-plane (φ=0°): slice at ky=0 ──────────────────────────────────────
Fx_ep = interp1(kx(:), fx(:,ky0), k0*sin_th, 'linear', 0);
Fy_ep = interp1(kx(:), fy(:,ky0), k0*sin_th, 'linear', 0);

E_theta_ep = cos_th .* Fx_ep;         % cosφ=1, sinφ=0
E_phi_ep   = Fy_ep;
E_tot_ep   = sqrt(abs(E_theta_ep).^2 + abs(E_phi_ep).^2);
E_ep_dB    = max(20*log10(E_tot_ep ./ max(E_tot_ep) + 1e-10), -40);

% ── H-plane (φ=90°): slice at kx=0 ─────────────────────────────────────
Fx_hp = interp1(ky(:), fx(kx0,:).', k0*sin_th, 'linear', 0);
Fy_hp = interp1(ky(:), fy(kx0,:).', k0*sin_th, 'linear', 0);

E_theta_hp = cos_th .* Fy_hp;         % cosφ=0, sinφ=1
E_phi_hp   = -Fx_hp;
E_tot_hp   = sqrt(abs(E_theta_hp).^2 + abs(E_phi_hp).^2);
E_hp_dB    = max(20*log10(E_tot_hp ./ max(E_tot_hp) + 1e-10), -40);

% ── Beamwidths ──────────────────────────────────────────────────────────
bw3_ep  = bw_calc(theta_deg, E_ep_dB, -3);
bw10_ep = bw_calc(theta_deg, E_ep_dB, -10);
bw3_hp  = bw_calc(theta_deg, E_hp_dB, -3);
bw10_hp = bw_calc(theta_deg, E_hp_dB, -10);

fprintf('\nE-plane (phi=0  ):  3-dB = %.1f deg  |  10-dB = %.1f deg\n', bw3_ep,  bw10_ep);
fprintf('H-plane (phi=90 ):  3-dB = %.1f deg  |  10-dB = %.1f deg\n',  bw3_hp,  bw10_hp);

%% ═══════════════════════════════════════════════════════════════════════════
%%  FIGURE 1 — Back-Projection  (2D colour maps, lecture note pg.41 style)
%%  Shows Ex and Ey: amplitude [dB] and phase [deg] → 4 subplots
%% ═══════════════════════════════════════════════════════════════════════════

% Back-projected complex fields (before abs), cropped to (M×N)
[Ex_c, Ey_c, ~, x_ax, y_ax] = ...
    BackProjection_PlanarNearField(E_meas_x, E_meas_y, z_0, f, dx, dy);
Ex_c = fliplr(Ex_c');   % (N×M) complex, viewed from front
Ey_c = fliplr(Ey_c');
Ex_norm = abs(Ex_c) ./ max(abs(Ex_c(:)));
Ey_norm = abs(Ey_c) ./ max(abs(Ey_c(:)));
Ex_dB   = 20*log10(Ex_norm + 1e-10);
Ey_dB   = 20*log10(Ey_norm + 1e-10);

fig1 = figure('Name','Back-Projection','Position',[50 30 1100 600]);
xm = x_ax*1e3;   ym = y_ax*1e3;   % [mm]

subplot(2,2,1);
imagesc(xm, ym, Ex_dB, [-40 0]); set(gca,'YDir','normal');
colormap(jet); colorbar; axis tight;
xlabel('x  (m)','FontSize',11); ylabel('y  (m)','FontSize',11);
title('E_x(abs)  [dB]','FontSize',12,'FontWeight','bold');
set(gca,'XTickLabel', arrayfun(@(v)sprintf('%.1f',v/1e3), ...
    get(gca,'XTick'),'UniformOutput',false));
set(gca,'YTickLabel', arrayfun(@(v)sprintf('%.1f',v/1e3), ...
    get(gca,'YTick'),'UniformOutput',false));

subplot(2,2,2);
imagesc(xm, ym, Ey_dB, [-40 0]); set(gca,'YDir','normal');
colormap(jet); colorbar; axis tight;
xlabel('x  (m)','FontSize',11); ylabel('y  (m)','FontSize',11);
title('E_y(abs)  [dB]','FontSize',12,'FontWeight','bold');
set(gca,'XTickLabel', arrayfun(@(v)sprintf('%.1f',v/1e3), ...
    get(gca,'XTick'),'UniformOutput',false));
set(gca,'YTickLabel', arrayfun(@(v)sprintf('%.1f',v/1e3), ...
    get(gca,'YTick'),'UniformOutput',false));

subplot(2,2,3);
imagesc(xm, ym, rad2deg(angle(Ex_c)), [-180 180]);
set(gca,'YDir','normal'); colormap(hsv); colorbar; axis tight;
xlabel('x  (m)','FontSize',11); ylabel('y  (m)','FontSize',11);
title('E_x(phase)  [deg]','FontSize',12,'FontWeight','bold');
set(gca,'XTickLabel', arrayfun(@(v)sprintf('%.1f',v/1e3), ...
    get(gca,'XTick'),'UniformOutput',false));
set(gca,'YTickLabel', arrayfun(@(v)sprintf('%.1f',v/1e3), ...
    get(gca,'YTick'),'UniformOutput',false));

subplot(2,2,4);
imagesc(xm, ym, rad2deg(angle(Ey_c)), [-180 180]);
set(gca,'YDir','normal'); colormap(hsv); colorbar; axis tight;
xlabel('x  (m)','FontSize',11); ylabel('y  (m)','FontSize',11);
title('E_y(phase)  [deg]','FontSize',12,'FontWeight','bold');
set(gca,'XTickLabel', arrayfun(@(v)sprintf('%.1f',v/1e3), ...
    get(gca,'XTick'),'UniformOutput',false));
set(gca,'YTickLabel', arrayfun(@(v)sprintf('%.1f',v/1e3), ...
    get(gca,'YTick'),'UniformOutput',false));

sgtitle({'SGH Back-Projection  @  5.3 GHz', ...
         sprintf('z_0 = %.0f mm  (Inverse 2D-FFT)', z_0*1e3)}, 'FontSize',13);

print(fig1, '-dpng', '-r150', 'project2_backproj.png');
fprintf('\n[Saved]  project2_backproj.png\n');

%% ═══════════════════════════════════════════════════════════════════════════
%%  FIGURE 2 — 1-D Far-Field Pattern  (Lecture Note 11.1, pg.40 style)
%%
%%  LEFT : H-plane (φ=90°)
%%    co-pol  = E_φ = −fx(0, k₀sinθ)              [main component]
%%    cross-pol = E_θ = cosθ·fy(0, k₀sinθ)
%%
%%  RIGHT: E-plane (φ=0°)
%%    co-pol  = E_θ = cosθ·fx(k₀sinθ, 0)          [main component]
%%    cross-pol = E_φ = fy(k₀sinθ, 0)
%% ═══════════════════════════════════════════════════════════════════════════
fig2 = figure('Name','Far-Field 1-D Pattern','Position',[100 80 1200 520]);

clip_dB = -80;

% ── H-plane co/cross-pol ────────────────────────────────────────────────
E_hp_copol   = to_dB_norm(abs(E_phi_hp),   clip_dB);
E_hp_xpol    = to_dB_norm(abs(E_theta_hp), clip_dB);

% ── E-plane co/cross-pol ────────────────────────────────────────────────
E_ep_copol   = to_dB_norm(abs(E_theta_ep), clip_dB);
E_ep_xpol    = to_dB_norm(abs(E_phi_ep),   clip_dB);

bw3_ep_co  = bw_calc(theta_deg, E_ep_copol, -3);
bw3_hp_co  = bw_calc(theta_deg, E_hp_copol, -3);

%% ── Subplot (left): H-plane ─────────────────────────────────────────────
subplot(1,2,1);
plot(theta_deg, E_hp_copol, 'r-o', 'LineWidth',1.5, ...
     'MarkerSize',3, 'MarkerIndices',1:50:numel(theta_deg));  hold on;
plot(theta_deg, E_hp_xpol,  'g-s', 'LineWidth',1.5, ...
     'MarkerSize',3, 'MarkerIndices',26:50:numel(theta_deg));
plot([-90 90],[-3  -3],  'k:', 'LineWidth',1.0);
plot([-90 90],[-10 -10], 'k--','LineWidth',1.0);
hold off;
xlim([-90 90]);  ylim([clip_dB 2]);
set(gca,'XTick',-90:20:90,'FontSize',11,'YDir','normal');
xlabel('\theta  [deg]','FontSize',13);
ylabel('normalized pattern  [dB]  (H-plane)','FontSize',12);
title({sprintf('E_{\\phi}^{far}(r,\\theta,\\phi) \\propto ik\\frac{e^{ikr}}{2\\pi r}\\cos\\theta\\left[-f_x\\sin\\phi+f_y\\cos\\phi\\right]'), ...
       sprintf('H-plane (\\phi=90°)   f = %.1f GHz   BW_3=%.1f°', f/1e9, bw3_hp_co)}, ...
      'FontSize',10);
legend('FFT (co-pol.)', 'FFT (cross-pol.)', '-3 dB', '-10 dB', ...
       'Location','best','FontSize',10);
grid on;

%% ── Subplot (right): E-plane ────────────────────────────────────────────
subplot(1,2,2);
plot(theta_deg, E_ep_copol, 'r-o', 'LineWidth',1.5, ...
     'MarkerSize',3, 'MarkerIndices',1:50:numel(theta_deg));  hold on;
plot(theta_deg, E_ep_xpol,  'g-s', 'LineWidth',1.5, ...
     'MarkerSize',3, 'MarkerIndices',26:50:numel(theta_deg));
plot([-90 90],[-3  -3],  'k:', 'LineWidth',1.0);
plot([-90 90],[-10 -10], 'k--','LineWidth',1.0);
hold off;
xlim([-90 90]);  ylim([clip_dB 2]);
set(gca,'XTick',-90:20:90,'FontSize',11,'YDir','normal');
xlabel('\theta  [deg]','FontSize',13);
ylabel('normalized pattern  [dB]  (E-plane)','FontSize',12);
title({sprintf('E_{\\theta}^{far}(r,\\theta,\\phi) \\propto ik\\frac{e^{ikr}}{2\\pi r}\\left[f_x\\cos\\phi+f_y\\sin\\phi\\right]'), ...
       sprintf('E-plane (\\phi=0°)   f = %.1f GHz   BW_3=%.1f°', f/1e9, bw3_ep_co)}, ...
      'FontSize',10);
legend('FFT (co-pol.)', 'FFT (cross-pol.)', '-3 dB', '-10 dB', ...
       'Location','best','FontSize',10);
grid on;

print(fig2, '-dpng', '-r150', 'project2_farfield.png');
fprintf('[Saved]  project2_farfield.png\n');

%% ═══════════════════════════════════════════════════════════════════════════
%%  FIGURE 3 — Plane Wave Spectrum  (auxiliary)
%% ═══════════════════════════════════════════════════════════════════════════
fig3 = figure('Name','Plane Wave Spectrum','Position',[200 200 580 490]);

kvis = 1.5 * k0;
kx_s = kx(abs(kx) <= kvis);
ky_s = ky(abs(ky) <= kvis);
pws_dB = 20*log10(abs(fx(abs(kx)<=kvis, abs(ky)<=kvis)) ./ ...
                   max(max(abs(fx))) + 1e-10);

imagesc(ky_s, kx_s, pws_dB, [-40 0]);
set(gca,'YDir','normal');
colormap(jet); colorbar;
hold on;
th_r = linspace(0, 2*pi, 361);
plot(k0*sin(th_r), k0*cos(th_r), 'w--', 'LineWidth', 1.8);
hold off;
xlabel('k_y  (rad/m)', 'FontSize',12);
ylabel('k_x  (rad/m)', 'FontSize',12);
title({'Plane Wave Spectrum  |f_x(k_x,k_y)|', ...
       'white dashed = visible region boundary  (k_\rho = k_0)'}, 'FontSize',11);
axis equal tight;

print(fig3, '-dpng', '-r150', 'project2_pws.png');
fprintf('[Saved]  project2_pws.png\n');

%% ═══ Local Functions ══════════════════════════════════════════════════════

function bw = bw_calc(theta, pat_dB, level)
    idx = find(pat_dB >= level);
    if isempty(idx), bw = NaN; return; end
    bw = theta(idx(end)) - theta(idx(1));
end

function y = to_dB_norm(E, clip)
    y = max(20*log10(abs(E)./max(abs(E(:)))+1e-12), clip);
end
