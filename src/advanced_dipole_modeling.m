%% advanced_dipole_modeling.m
%  Advanced Project: Infinitesimal Dipole Modeling of Horn Antenna
%  Course  : Defense Antenna Engineering
%  Student : 202550592 오선지
%  AUT     : Standard Gain Horn (SGH) @ 5.3 GHz
%
%  Theory  : Lecture Note 14.1 / Balanis §3.6
%  Method  : Near-field → Hertzian dipole array (inverse problem)
%             → reconstruct far-field pattern & aperture image
%
%  Outputs:
%    advanced_aperture.png       — NF fit & aperture images
%    advanced_dipole_results.png — far-field comparison (Dipole vs PWS)

clear; clc; close all;

%% ─── Path & Data ────────────────────────────────────────────────────────
script_dir = fileparts(mfilename('fullpath'));
data_dir   = fullfile(script_dir, 'pr2_files');
addpath(data_dir);

load(fullfile(data_dir, 'SGH_FullPlanarScan.mat'));

%% ─── Parameters ─────────────────────────────────────────────────────────
c     = 299792458;
f     = 5.3e9;
lam   = c / f;
k0    = 2*pi / lam;
omega = 2*pi * f;
mu0   = 4*pi * 1e-7;
z_0   = 0.088;
dx    = 0.022;   dy = 0.022;

E_meas_x = comp_x_measurement(:,:,3);   % (68×36)
E_meas_y = comp_y_measurement(:,:,3);
[M, N]   = size(E_meas_x);

% Centred measurement grid
xm = ((0:M-1) - M/2) * dx;   % (1×M)
ym = ((0:N-1) - N/2) * dy;   % (1×N)
[XM_g, YM_g] = ndgrid(xm, ym);   % (M×N), same as Python meshgrid(indexing='ij')

%% ═══ Dipole Grid ══════════════════════════════════════════════════════════
%  Coarser than measurement (stride=2) → avoids over-fitting
%  34 × 18 = 612 dipoles

S   = 2;
Ndx = M / S;            % 34
Ndy = N / S;            % 18
Nd  = Ndx * Ndy;        % 612

xd  = ((0:Ndx-1) - Ndx/2) * S * dx;   % (1×Ndx)
yd  = ((0:Ndy-1) - Ndy/2) * S * dy;   % (1×Ndy)
[XD, YD] = ndgrid(xd, yd);            % (Ndx×Ndy)
xd_f = XD(:);                          % (Nd×1)
yd_f = YD(:);

fprintf('Scan   : %d x %d = %d points\n', M, N, M*N);
fprintf('Dipoles: %d x %d = %d  (stride=%d)\n', Ndx, Ndy, Nd, S);
fprintf('Unknowns  = %d  |  Equations = %d\n', 2*Nd, 2*M*N);

%% ═══ Green's Function Matrix  G  (2Nm × 2Nd) ═════════════════════════════
%
%  x-dipole at r'=(x',y',0), observation (x,y,z0)  [Balanis §3.6]:
%    E_x = C·g·[(1−Rx²) + (3Rx²−1)·(j/kR + 1/(kR)²)]
%    E_y = C·g·Rx·Ry·[−1 + 3j/kR + 3/(kR)²]
%  y-dipole: swap x↔y  (by symmetry Gyx = Gxy)
%    C = −jωμ₀,   g = exp(−jkR)/(4πR)

fprintf('\nBuilding Green function matrix G ...\n');
tic;
G = build_G(xm, ym, xd_f, yd_f, z_0, k0, omega, mu0, XM_g, YM_g);
fprintf('G : %d x %d  (%.0f MB)  [%.1f s]\n', ...
        size(G,1), size(G,2), numel(G)*16/1e6, toc);

%% ═══ Tikhonov Regularised Inverse  p = (G^H G + λI)^{-1} G^H b ══════════
Nm  = M * N;
b   = [E_meas_x(:); E_meas_y(:)];   % (2Nm×1)

GHG   = G' * G;                      % conjugate transpose
GHb   = G' * b;
alpha = 1e-2;
lam_r = alpha * max(real(diag(GHG)));
fprintf('Regularisation  lambda = %.3e  (alpha=%.0e)\n', lam_r, alpha);

fprintf('Solving ...\n');
tic;
p_opt = (GHG + lam_r * eye(2*Nd)) \ GHb;
fprintf('Done in %.1f s\n', toc);

px_f = p_opt(1:Nd);          % x-dipole moments  (Nd×1)
py_f = p_opt(Nd+1:end);      % y-dipole moments

b_pred = G * p_opt;
resid  = norm(b - b_pred) / norm(b);
fprintf('Relative residual : %.1f%%\n\n', resid*100);

%% ═══ Far-Field Pattern — Dipole Model ═════════════════════════════════════
%
%  E_θ(θ,φ) = Σⱼ exp(jk·r̂·rⱼ) · [pxⱼ cosθ cosφ + pyⱼ cosθ sinφ]
%  E_φ(θ,φ) = Σⱼ exp(jk·r̂·rⱼ) · [−pxⱼ sinφ + pyⱼ cosφ]

theta_deg = linspace(-90, 90, 1801);
theta_rad = deg2rad(theta_deg);
sin_th    = sin(theta_rad);    % (1×Nθ) row
cos_th    = cos(theta_rad);

E_dip_ep = ff_pattern(0,  xd_f, yd_f, px_f, py_f, k0, sin_th, cos_th);
E_dip_hp = ff_pattern(90, xd_f, yd_f, px_f, py_f, k0, sin_th, cos_th);
E_dip_ep_dB = to_dB(E_dip_ep);
E_dip_hp_dB = to_dB(E_dip_hp);

%% ═══ Far-Field Pattern — PWS Reference (Project 2) ════════════════════════
MI = 10*M;  NI = 10*N;
fx = ifftshift(ifft2(E_meas_x, MI, NI));
fy = ifftshift(ifft2(E_meas_y, MI, NI));
m  = (-MI/2 : MI/2-1);   kx = 2*pi*m/(MI*dx);
n  = (-NI/2 : NI/2-1);   ky = 2*pi*n/(NI*dy);

ky0 = NI/2+1;  kx0 = MI/2+1;
Fx_ep = interp1(kx(:), fx(:,ky0),   k0*sin_th, 'linear', 0);
Fy_ep = interp1(kx(:), fy(:,ky0),   k0*sin_th, 'linear', 0);
Fx_hp = interp1(ky(:), fx(kx0,:).', k0*sin_th, 'linear', 0);
Fy_hp = interp1(ky(:), fy(kx0,:).', k0*sin_th, 'linear', 0);

E_pws_ep = sqrt(abs(cos_th.*Fx_ep).^2 + abs(Fy_ep).^2);
E_pws_hp = sqrt(abs(cos_th.*Fy_hp).^2 + abs(Fx_hp).^2);
E_pws_ep_dB = to_dB(E_pws_ep);
E_pws_hp_dB = to_dB(E_pws_hp);

%% ── Beamwidths ────────────────────────────────────────────────────────────
bw_dip_ep  = bw_calc(theta_deg, E_dip_ep_dB, -3);
bw_dip_hp  = bw_calc(theta_deg, E_dip_hp_dB, -3);
bw_pws_ep  = bw_calc(theta_deg, E_pws_ep_dB, -3);
bw_pws_hp  = bw_calc(theta_deg, E_pws_hp_dB, -3);

fprintf('Dipole model:  E-plane BW3=%.1f°  H-plane BW3=%.1f°\n', bw_dip_ep, bw_dip_hp);
fprintf('PWS method:    E-plane BW3=%.1f°  H-plane BW3=%.1f°\n', bw_pws_ep, bw_pws_hp);

%% ═══ Aperture Images ══════════════════════════════════════════════════════
Ex_pred = reshape(G(1:Nm,:)    * p_opt, M, N);
Ey_pred = reshape(G(Nm+1:end,:)* p_opt, M, N);
U_pred  = abs(Ex_pred).^2 + abs(Ey_pred).^2;
U_meas  = abs(E_meas_x).^2  + abs(E_meas_y).^2;
U_dip   = reshape(abs(px_f).^2 + abs(py_f).^2, Ndx, Ndy);

% Back-projection (PWS)
[KY_G, KX_G] = meshgrid(ky, kx);
KZ_G      = sqrt(k0^2 - KX_G.^2 - KY_G.^2 + 0j);
prop_mask = (KX_G.^2 + KY_G.^2) < k0^2;
fx_bp = zeros(MI,NI);  fy_bp = zeros(MI,NI);
fx_bp(prop_mask) = fx(prop_mask) .* exp(1j*KZ_G(prop_mask)*z_0);
fy_bp(prop_mask) = fy(prop_mask) .* exp(1j*KZ_G(prop_mask)*z_0);
Ex_bp = abs(fft2(ifftshift(fx_bp)));  Ex_bp = Ex_bp(1:M,1:N);
Ey_bp = abs(fft2(ifftshift(fy_bp)));  Ey_bp = Ey_bp(1:M,1:N);
U_bp  = Ex_bp.^2 + Ey_bp.^2;

xm_mm = xm*1e3;   ym_mm = ym*1e3;
xd_mm = xd*1e3;   yd_mm = yd*1e3;

%% ═══════════════════════════════════════════════════════════════════════════
%%  FIGURE 1 — NF Fit & Aperture Images
%% ═══════════════════════════════════════════════════════════════════════════
fig1 = figure('Name','Dipole Modeling — Aperture','Position',[50 50 1300 470]);

subplot(1,4,1);
imagesc(ym_mm, xm_mm, norm_dB(U_meas));
set(gca,'YDir','normal'); colormap(jet); colorbar; clim([-40 0]);
xlabel('y (mm)'); ylabel('x (mm)');
title('Measured NF  |E_x|^2+|E_y|^2');

subplot(1,4,2);
imagesc(ym_mm, xm_mm, norm_dB(U_pred));
set(gca,'YDir','normal'); colormap(jet); colorbar; clim([-40 0]);
xlabel('y (mm)'); ylabel('x (mm)');
title(sprintf('Dipole-Predicted NF\n(residual %.1f%%)', resid*100));

subplot(1,4,3);
imagesc(yd_mm, xd_mm, norm_dB(U_dip));
set(gca,'YDir','normal'); colormap(jet); colorbar; clim([-40 0]);
xlabel('y (mm)'); ylabel('x (mm)');
title(sprintf('Dipole Aperture Image\n%dx%d=%d dipoles', Ndx, Ndy, Nd));

subplot(1,4,4);
imagesc(ym_mm, xm_mm, norm_dB(U_bp));
set(gca,'YDir','normal'); colormap(jet); colorbar; clim([-40 0]);
xlabel('y (mm)'); ylabel('x (mm)');
title('Back-Projection (PWS, Proj.2)');

sgtitle(sprintf('SGH @ %.1f GHz  —  Advanced: Dipole Modeling vs PWS', f/1e9), ...
        'FontSize',13);

print(fig1, '-dpng', '-r150', 'advanced_aperture.png');
fprintf('\n[Saved]  advanced_aperture.png\n');

%% ═══════════════════════════════════════════════════════════════════════════
%%  FIGURE 2 — Far-Field Pattern Comparison (Dipole vs PWS)
%% ═══════════════════════════════════════════════════════════════════════════
fig2 = figure('Name','Far-Field — Dipole vs PWS','Position',[100 100 1050 460]);

%% ── (a) Cartesian ──────────────────────────────────────────────────────
subplot(1,2,1);
plot(theta_deg, E_dip_ep_dB, 'b-',  'LineWidth', 2.2);  hold on;
plot(theta_deg, E_dip_hp_dB, 'r-',  'LineWidth', 2.2);
plot(theta_deg, E_pws_ep_dB, 'b--', 'LineWidth', 1.5);
plot(theta_deg, E_pws_hp_dB, 'r--', 'LineWidth', 1.5);
plot([-90 90],[-3  -3],  'k:',  'LineWidth', 1.0);
plot([-90 90],[-10 -10], 'k-.', 'LineWidth', 1.0);
hold off;
xlim([-90 90]);  ylim([-40 2]);
set(gca,'XTick',-90:15:90,'FontSize',11);
xlabel('Observation Angle  \theta  (deg)', 'FontSize',12);
ylabel('Normalized Pattern  (dB)',          'FontSize',12);
title(sprintf('Far-Field Pattern  (f = %.1f GHz)', f/1e9), 'FontSize',12);
legend(sprintf('Dipole  E-plane  BW_3=%.1f°', bw_dip_ep), ...
       sprintf('Dipole  H-plane  BW_3=%.1f°', bw_dip_hp), ...
       sprintf('PWS     E-plane  BW_3=%.1f°', bw_pws_ep), ...
       sprintf('PWS     H-plane  BW_3=%.1f°', bw_pws_hp), ...
       '-3 dB', '-10 dB', 'Location','best', 'FontSize',9);
grid on;

%% ── (b) Polar ──────────────────────────────────────────────────────────
subplot(1,2,2);
polarplot(theta_rad, max(E_dip_ep_dB+40,0), 'b-',  'LineWidth',2.2); hold on;
polarplot(theta_rad, max(E_dip_hp_dB+40,0), 'r-',  'LineWidth',2.2);
polarplot(theta_rad, max(E_pws_ep_dB+40,0), 'b--', 'LineWidth',1.5);
polarplot(theta_rad, max(E_pws_hp_dB+40,0), 'r--', 'LineWidth',1.5);
hold off;
ax = gca;
try; ax.ThetaZeroLocation = 'top';    catch; end
try; ax.ThetaDirection    = 'clockwise'; catch; end
try; ax.ThetaLim          = [-90 90]; catch; end
try; ax.RLim              = [0 42];   catch; end
try; ax.RTick             = [10 20 30 40]; catch; end
try; ax.RTickLabel        = {'-30 dB','-20 dB','-10 dB','0 dB'}; catch; end
try; ax.FontSize          = 10;       catch; end
legend('Dipole E','Dipole H','PWS E','PWS H', ...
       'Location','southoutside','Orientation','horizontal','FontSize',9);
title('Polar Pattern  (Dipole vs PWS)', 'FontSize',12);

print(fig2, '-dpng', '-r150', 'advanced_dipole_results.png');
fprintf('[Saved]  advanced_dipole_results.png\n');

fprintf('\n%s\n  ADVANCED — SUMMARY\n%s\n', repmat('=',1,50), repmat('=',1,50));
fprintf('  Dipole  E-plane 3-dB BW = %.1f deg\n', bw_dip_ep);
fprintf('  Dipole  H-plane 3-dB BW = %.1f deg\n', bw_dip_hp);
fprintf('  PWS     E-plane 3-dB BW = %.1f deg\n', bw_pws_ep);
fprintf('  PWS     H-plane 3-dB BW = %.1f deg\n', bw_pws_hp);
fprintf('  NF fit residual = %.1f%%\n', resid*100);
fprintf('%s\n', repmat('=',1,50));


%% ═══ Local Functions ══════════════════════════════════════════════════════

function G = build_G(xm, ym, xd_f, yd_f, z0, k0, omega, mu0, XM_g, YM_g)
% Green's function matrix  G  (2Nm × 2Nd)  [Balanis §3.6]
%
%  ndgrid convention (same as Python indexing='ij'):
%    XM_g, YM_g are (M×N) grids passed in to avoid re-computation
    xmf = XM_g(:);   ymf = YM_g(:);   % (Nm×1)

    % (Nm × Nd) — column vector minus row vector → matrix (implicit expansion)
    dX = xmf - xd_f(:).';    % (Nm, Nd)
    dY = ymf - yd_f(:).';
    R  = sqrt(dX.^2 + dY.^2 + z0^2);
    Rx = dX ./ R;
    Ry = dY ./ R;

    g  = exp(-1j*k0*R) ./ (4*pi*R);   % scalar Green's function
    kR = k0 * R;
    F1 = 1j./kR + 1./kR.^2;           % near-field term
    Fc = -1 + 3j./kR + 3./kR.^2;      % cross-coupling term
    C  = -1j * omega * mu0;

    Gxx = C*g.*((1-Rx.^2)     + (3*Rx.^2-1).*F1);   % x-dip → Ex
    Gxy = C*g.*Rx.*Ry.*Fc;                            % y-dip → Ex
    Gyx = Gxy;                                         % symmetry
    Gyy = C*g.*((1-Ry.^2)     + (3*Ry.^2-1).*F1);   % y-dip → Ey

    G = [Gxx, Gxy;
         Gyx, Gyy];      % (2Nm × 2Nd)
end


function E_tot = ff_pattern(phi_deg, xd_f, yd_f, px_f, py_f, k0, sin_th, cos_th)
% Far-field from dipole array  (array-factor phase summation)
%  xd_f, yd_f, px_f, py_f : (Nd×1)
%  sin_th, cos_th          : (1×Nθ)  →  implicit broadcasting → (Nd×Nθ)
    phi = deg2rad(phi_deg);
    cp  = cos(phi);   sp = sin(phi);

    phase = exp(1j*k0 * (xd_f.*sin_th.*cp + yd_f.*sin_th.*sp));  % (Nd×Nθ)

    E_th = sum(phase .* (px_f.*cos_th.*cp + py_f.*cos_th.*sp), 1);
    E_ph = sum(phase .* (-px_f.*sp         + py_f.*cp        ), 1);
    E_tot = sqrt(abs(E_th).^2 + abs(E_ph).^2);   % (1×Nθ)
end


function y = to_dB(E, clip)
    if nargin < 2, clip = -40; end
    y = max(20*log10(E ./ max(E(:)) + 1e-12), clip);
end

function y = norm_dB(A, clip)
    if nargin < 2, clip = -40; end
    y = max(10*log10(A ./ max(A(:)) + 1e-12), clip);
end

function bw = bw_calc(theta, pat_dB, level)
    idx = find(pat_dB >= level);
    if isempty(idx), bw = NaN; return; end
    bw = theta(idx(end)) - theta(idx(1));
end
