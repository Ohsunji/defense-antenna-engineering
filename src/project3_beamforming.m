%% Project 3: 안테나 배열 빔포밍 (Array Beamforming)
%  이론: 강의노트 7 (Antenna Array and Beamforming)
%  배열: 4×4 UPA (Uniform Planar Array), 소자 = isotropic
%
%  [핵심 행렬식 기반 공식]
%
%  스티어링 벡터 (N×1):
%    a(θ,φ) = exp(j·k₀·r · k̂)
%    r: 소자 위치 행렬 (N×3), k̂: 단위 방향 벡터 (3×1)
%
%  배열 인수 (Array Factor):
%    AF(θ,φ) = w^H · a(θ,φ)   ← 행렬 내적
%
%  스티어링 빔계수:
%    w = a(θ₀,φ₀)
%    → w^H·a(θ₀,φ₀) = a^H·a = Σ|a_n|² = N  (최대값, 완전 위상 정렬)
%
%  지향성:
%    D_max = 4π × P_max / ∫∫ P(θ,φ) sinθ dθ dφ

clear; clc; close all;

%% ===== 1. 배열 파라미터 =====
f      = 10e9;           % 주파수: 10 GHz
c      = 3e8;            % 광속 [m/s]
lambda = c / f;          % 파장 λ [m]
k0     = 2*pi / lambda;  % 파수 k₀ = 2π/λ [rad/m]
d      = lambda / 2;     % 소자 간격: λ/2 (grating lobe 방지 조건 d < λ)

Nx = 4;          % x방향 소자 수
Ny = 4;          % y방향 소자 수
N  = Nx * Ny;    % 총 소자 수 = 16

% 소자 위치 행렬 r (N×3): x-y 평면에 배치
% (mx·d, ny·d, 0), mx = 0..Nx-1, ny = 0..Ny-1
[mx, ny_idx] = meshgrid(0:Nx-1, 0:Ny-1);
r = [mx(:)*d,  ny_idx(:)*d,  zeros(N,1)];   % N×3

fprintf('===== 배열 파라미터 =====\n');
fprintf('  주파수:   %.0f GHz\n', f/1e9);
fprintf('  파장 λ:   %.4f m\n', lambda);
fprintf('  소자 간격: λ/2 = %.4f m\n', d);
fprintf('  배열 크기: %d×%d = %d 소자 (UPA)\n', Nx, Ny, N);
fprintf('  Grating lobe 조건: d=λ/2 < λ → 만족\n');

%% ===== 2. 각도 그리드 =====
nTheta = 181;   % 0° ~ 180°
nPhi   = 361;   % 0° ~ 360°

theta_vec = linspace(0,   pi,   nTheta);  % [rad]
phi_vec   = linspace(0, 2*pi,   nPhi  );  % [rad]

[Phi_g, Theta_g] = meshgrid(phi_vec, theta_vec);  % nTheta × nPhi

% 단위 방향 벡터 k̂ = [sin θ cos φ, sin θ sin φ, cos θ]
% K: 3 × (nTheta*nPhi) — 모든 방향의 k̂ 행렬
K = k0 * [sin(Theta_g(:)).*cos(Phi_g(:)), ...
           sin(Theta_g(:)).*sin(Phi_g(:)), ...
           cos(Theta_g(:))]';         % 3 × (nTheta*nPhi)

%% ===== 3. 스티어링 행렬 A 구성 (핵심 행렬 연산) =====
%
%  Phase_mat(n,m) = k₀ · r_n · k̂_m  (n: 소자, m: 방향)
%  A(n,m) = exp(j · Phase_mat(n,m))
%
%  A: N × (nTheta*nPhi)
%  각 열이 하나의 (θ,φ) 방향의 스티어링 벡터 a(θ,φ)

Phase_mat = r * K;                 % (N×3) × (3×nAngles) = N×nAngles
A = exp(1j * Phase_mat);           % N × (nTheta*nPhi)

fprintf('\n스티어링 행렬 A: %d × %d\n', size(A,1), size(A,2));

%% ===== 4. Case 1: 빔계수 모두 1 (Uniform Weighting) =====
%
%  w_uniform = [1, 1, ..., 1]^T  (N×1)
%  → 위상 정렬 없음, 기본 배열 인수

w1   = ones(N, 1);             % 빔계수: 모두 1
AF1  = w1' * A;                % (1×N) × (N×nAngles) = 1×nAngles
P1   = reshape(abs(AF1).^2, nTheta, nPhi);  % 전력 패턴

fprintf('\n[Case 1] 빔계수: w = ones(%d,1)\n', N);

%% ===== 5. Case 2: θ₀ 방향 스티어링 빔계수 추출 =====
%
%  목표: 특정 방향 (θ₀,φ₀)에서 AF가 최대 (위상 완전 정렬)
%
%  스티어링 벡터: a(θ₀,φ₀) = exp(j·k₀·r·k̂₀)
%  최적 빔계수:  w = a(θ₀,φ₀)
%
%  증명 (MATLAB: w'*a = Σ conj(w_n)*a_n):
%    AF(θ₀,φ₀) = w^H · a(θ₀,φ₀)
%               = a(θ₀,φ₀)^H · a(θ₀,φ₀)   [w = a이므로 w^H = a^H]
%               = Σ conj(a_n)·a_n = Σ|a_n|² = N  (최대)

theta0_deg = 30;   phi0_deg = 0;    % 스티어링 목표 방향
theta0 = deg2rad(theta0_deg);
phi0   = deg2rad(phi0_deg);

% 목표 방향 단위 벡터
k_hat0 = [sin(theta0)*cos(phi0); sin(theta0)*sin(phi0); cos(theta0)];

% 스티어링 벡터 a(θ₀,φ₀): N×1
a_steer = exp(1j * k0 * (r * k_hat0));   % N×1

% 빔계수 추출
%  MATLAB에서 w'*a = Σ conj(w_n)*a_n 이므로
%  w = a_steer 로 설정해야:
%    w^H * a(θ₀) = Σ conj(a_n)*a_n = Σ|a_n|² = N (최대) ✓
%  (w = conj(a_steer) 로 설정하면 Σ a_n² → φ=180° 방향으로 잘못 스티어링됨)
w2 = a_steer;   % N×1

% 검증: 스티어링 방향에서 AF = N (최대)
AF_check = w2' * a_steer;
fprintf('\n[Case 2] 스티어링 방향: θ₀=%d°, φ₀=%d°\n', theta0_deg, phi0_deg);
fprintf('  빔계수: w_n = exp(+j·k₀·r_n·k̂₀) = a_n(θ₀,φ₀)  (n=1..%d)\n', N);
fprintf('  검증: AF(θ₀,φ₀) = %.4f  (이론값 N=%d) ✓\n', abs(AF_check), N);

% 배열 인수 및 전력 패턴
AF2 = w2' * A;
P2  = reshape(abs(AF2).^2, nTheta, nPhi);

%% ===== 6. 지향성(Directivity) 계산 =====
%
%  D(θ₀,φ₀) = 4π × U_max / P_rad
%  P_rad = ∫∫ U(θ,φ) sinθ dθ dφ ≈ ΣΣ P(θ,φ)·sinθ·Δθ·Δφ

dTheta = theta_vec(2) - theta_vec(1);  % [rad]
dPhi   = phi_vec(2)   - phi_vec(1);    % [rad]
sin_th = sin(theta_vec)';              % nTheta×1

P_rad1 = sum(sum(P1 .* sin_th)) * dTheta * dPhi;
D1_lin = 4*pi * max(P1(:)) / P_rad1;
D1_dBi = 10*log10(D1_lin);

P_rad2 = sum(sum(P2 .* sin_th)) * dTheta * dPhi;
D2_lin = 4*pi * max(P2(:)) / P_rad2;
D2_dBi = 10*log10(D2_lin);

fprintf('\n===== 지향성(Directivity) =====\n');
fprintf('  [Case 1] Uniform:  D = %.2f (선형) = %.2f dBi\n', D1_lin, D1_dBi);
fprintf('  [Case 2] Steered:  D = %.2f (선형) = %.2f dBi\n', D2_lin, D2_dBi);
fprintf('  (이론 참고: 최대 D ≈ N = %d = %.1f dBi)\n', N, 10*log10(N));

%% ===== 7. 1D 패턴: E-Plane (φ=0°) & H-Plane (φ=90°) =====
phi0_idx  = 1;
phi90_idx = find(phi_vec >= pi/2, 1);

theta_deg = rad2deg(theta_vec);
dB_norm   = @(P) 10*log10(P ./ max(P(:)) + 1e-10);  % 정규화 dB (각 컷의 local max 기준)

figure('Name','Project3 - 1D Patterns','Position',[50 50 1200 520]);

% --- E-Plane (φ=0°) ---
subplot(1,2,1);
plot(theta_deg, dB_norm(P1(:,phi0_idx)), 'b-',  'LineWidth', 2.0); hold on;
plot(theta_deg, dB_norm(P2(:,phi0_idx)), 'r-',  'LineWidth', 2.0);
xline(theta0_deg, 'k--', sprintf('\\theta_0=%d°', theta0_deg), ...
      'LineWidth', 1.2, 'LabelVerticalAlignment', 'bottom');
grid on;
xlabel('Theta (deg)', 'FontSize', 11);
ylabel('Normalized Gain (dB)', 'FontSize', 11);
title('E-Plane  (φ = 0°)', 'FontSize', 12);
legend('Uniform (w=1)', sprintf('Steered (\\theta_0=%d°)', theta0_deg), ...
       'Location', 'best', 'FontSize', 10);
xlim([0 180]); ylim([-40 5]);

% --- H-Plane (φ=90°) ---
subplot(1,2,2);
plot(theta_deg, dB_norm(P1(:,phi90_idx)), 'b-',  'LineWidth', 2.0); hold on;
plot(theta_deg, dB_norm(P2(:,phi90_idx)), 'r-',  'LineWidth', 2.0);
grid on;
xlabel('Theta (deg)', 'FontSize', 11);
ylabel('Normalized Gain (dB)', 'FontSize', 11);
title('H-Plane  (φ = 90°)', 'FontSize', 12);
legend('Uniform (w=1)', sprintf('Steered (\\theta_0=%d°)', theta0_deg), ...
       'Location', 'best', 'FontSize', 10);
xlim([0 180]); ylim([-40 5]);

sgtitle(sprintf('%d×%d UPA, d=\\lambda/2, f=%.0f GHz  |  D_{uni}=%.1f dBi, D_{steer}=%.1f dBi', ...
    Nx, Ny, f/1e9, D1_dBi, D2_dBi), 'FontSize', 12);

%% ===== 7-2. 1D 패턴: -180°~180° (E-Plane & H-Plane) =====
%
%  구성:
%    angle < 0 → phi=180°(E) / phi=270°(H) 컷, theta 180°→0°
%    angle ≥ 0 → phi=0°(E)  / phi=90°(H)  컷, theta 0°→180°

phi180_idx = find(phi_vec >= pi,     1);   % φ=180°
phi270_idx = find(phi_vec >= 3*pi/2, 1);   % φ=270°

angle_axis = linspace(-180, 180, 2*nTheta-1);  % -180:1:180 (361점)

% ── E-Plane 합성 ──
% 왼쪽(-180→0): phi=180°, theta 180→0 = flipud
E1_left = flipud(dB_norm(P1(:, phi180_idx)));
E2_left = flipud(dB_norm(P2(:, phi180_idx)));
% 오른쪽(1→180): phi=0°, theta 1→180 (theta=0 제외로 중복 방지)
E1_right = dB_norm(P1(2:end, phi0_idx));
E2_right = dB_norm(P2(2:end, phi0_idx));

E1_full = [E1_left; E1_right];
E2_full = [E2_left; E2_right];

% ── H-Plane 합성 ──
H1_left = flipud(dB_norm(P1(:, phi270_idx)));
H2_left = flipud(dB_norm(P2(:, phi270_idx)));
H1_right = dB_norm(P1(2:end, phi90_idx));
H2_right = dB_norm(P2(2:end, phi90_idx));

H1_full = [H1_left; H1_right];
H2_full = [H2_left; H2_right];

figure('Name','Project3 - 1D Patterns (-180~180)','Position',[50 50 1200 520]);

% --- E-Plane ---
subplot(1,2,1);
plot(angle_axis, E1_full, 'b-', 'LineWidth', 2.0); hold on;
plot(angle_axis, E2_full, 'r-', 'LineWidth', 2.0);
xline(theta0_deg,  'k--', sprintf('\\theta_0=%d°', theta0_deg),  ...
      'LineWidth', 1.2, 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
xline(-theta0_deg, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
xline(0, 'k:', 'LineWidth', 1.0, 'HandleVisibility', 'off');
grid on;
xlabel('Angle (deg)', 'FontSize', 11);
ylabel('Normalized Gain (dB)', 'FontSize', 11);
title('E-Plane  (φ = 0°/180°)', 'FontSize', 12);
legend('Uniform (w=1)', sprintf('Steered (\\theta_0=%d°)', theta0_deg), ...
       'Location', 'best', 'FontSize', 10);
xlim([-180 180]); xticks(-180:30:180); ylim([-40 5]);

% --- H-Plane ---
subplot(1,2,2);
plot(angle_axis, H1_full, 'b-', 'LineWidth', 2.0); hold on;
plot(angle_axis, H2_full, 'r-', 'LineWidth', 2.0);
xline(0, 'k:', 'LineWidth', 1.0, 'HandleVisibility', 'off');
grid on;
xlabel('Angle (deg)', 'FontSize', 11);
ylabel('Normalized Gain (dB)', 'FontSize', 11);
title('H-Plane  (φ = 90°/270°)', 'FontSize', 12);
legend('Uniform (w=1)', sprintf('Steered (\\theta_0=%d°)', theta0_deg), ...
       'Location', 'best', 'FontSize', 10);
xlim([-180 180]); xticks(-180:30:180); ylim([-40 5]);

sgtitle(sprintf('%d×%d UPA, d=\\lambda/2, f=%.0f GHz  |  D_{uni}=%.1f dBi, D_{steer}=%.1f dBi', ...
    Nx, Ny, f/1e9, D1_dBi, D2_dBi), 'FontSize', 12);

%% ===== 8. 3D 패턴 시각화 =====
%  구면좌표 → 직교좌표: x=r·sinθ·cosφ, y=r·sinθ·sinφ, z=r·cosθ
%  색상 = 정규화 이득 [dB]

titles_3d = {'Case 1: Uniform (w = 1)', ...
             sprintf('Case 2: Steered (\\theta_0=%d°, \\phi_0=%d°)', theta0_deg, phi0_deg)};
P_list    = {P1, P2};
D_list    = [D1_dBi, D2_dBi];

for ci = 1:2
    P     = P_list{ci};
    R_n   = P / max(P(:));                   % 정규화 (0~1)
    X = R_n .* sin(Theta_g) .* cos(Phi_g);
    Y = R_n .* sin(Theta_g) .* sin(Phi_g);
    Z = R_n .* cos(Theta_g);
    C = dB_norm(P);                           % 색상: dB값

    figure('Name', sprintf('Project3 - 3D (%s)', titles_3d{ci}), ...
           'Position', [50+480*(ci-1), 620, 700, 600]);
    surf(X, Y, Z, C, 'EdgeColor', 'none');
    colorbar; colormap jet; clim([-30 0]);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(sprintf('%s\nD_{max} = %.2f dBi', titles_3d{ci}, D_list(ci)), ...
          'FontSize', 11);
    axis equal; grid on; view(-45, 30);
end

%% ===== 9. 빔계수 시각화 (크기 & 위상) =====
figure('Name','Project3 - Beamweights','Position',[50 50 1000 400]);

w2_mag   = abs(w2);
w2_phase = rad2deg(angle(w2));

% 소자 인덱스 (2D 배치)
elem_idx = reshape(1:N, Ny, Nx);

subplot(1,2,1);
imagesc(reshape(w2_mag, Ny, Nx));
colorbar; colormap hot;
xlabel('X 소자 인덱스'); ylabel('Y 소자 인덱스');
title('빔계수 크기 |w_n|  (모두 1.0)');
axis equal tight;

subplot(1,2,2);
imagesc(reshape(w2_phase, Ny, Nx));
colorbar; colormap hsv;
xlabel('X 소자 인덱스'); ylabel('Y 소자 인덱스');
title(sprintf('빔계수 위상 ∠w_n [deg]  (θ₀=%d°, φ₀=%d°)', theta0_deg, phi0_deg));
axis equal tight;

sgtitle(sprintf('스티어링 빔계수: w_n = exp(-j·k₀·r_n·k̂₀)'), 'FontSize', 12);

%% ===== 결과 요약 =====
fprintf('\n========== Project 3 결과 요약 ==========\n');
fprintf('  배열:          %d×%d UPA, d=λ/2, f=%.0f GHz\n', Nx, Ny, f/1e9);
fprintf('  [Uniform]  D = %.2f dBi\n', D1_dBi);
fprintf('  [Steered]  D = %.2f dBi  @ θ₀=%d°\n', D2_dBi, theta0_deg);
fprintf('  이론값: 최대 D ≈ %d = %.1f dBi  (N개 isotropic 소자)\n', N, 10*log10(N));
fprintf('==========================================\n');
