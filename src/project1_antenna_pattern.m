%% Project 1: 안테나 패턴 3D 출력 및 지향성(Directivity) 계산
%  이론: 강의노트 5 (Antenna Measurement - Far Field)
%  데이터: pr1_far_field_mat.mat
%
%  데이터 구조 (7381행 × 6열):
%    Col 1: Theta [deg], 범위 0~180
%    Col 2: Phi   [deg], 범위 -180~180
%    Col 3: Vert  Gain  (E_theta 성분) [dBi]
%    Col 4: Vert  Phase (E_theta 성분) [deg]
%    Col 5: Horiz Gain  (E_phi   성분) [dBi]
%    Col 6: Horiz Phase (E_phi   성분) [deg]
%  정렬 순서: phi 외부 루프(-180→180, 3°간격), theta 내부 루프(180→0, 3°간격)

clear; clc; close all;

%% ===== 1. 데이터 로드 =====
load('pr1_far_field_mat.mat');  % 변수명: far_field_mat

nTheta = 61;   % theta: 180→0 (3°간격, 61개)
nPhi   = 121;  % phi: -180→180 (3°간격, 121개)

% 1D 벡터에서 2D 행렬로 reshape
% 결과: 행=theta 방향(180→0), 열=phi 방향(-180→180)
theta_mat    = reshape(far_field_mat(:,1), nTheta, nPhi);  % [deg]
phi_mat      = reshape(far_field_mat(:,2), nTheta, nPhi);  % [deg]
Gv_dBi_mat   = reshape(far_field_mat(:,3), nTheta, nPhi); % E_theta 이득 [dBi]
Pv_deg_mat   = reshape(far_field_mat(:,4), nTheta, nPhi); % E_theta 위상 [deg]
Gh_dBi_mat   = reshape(far_field_mat(:,5), nTheta, nPhi); % E_phi 이득 [dBi]
Ph_deg_mat   = reshape(far_field_mat(:,6), nTheta, nPhi); % E_phi 위상 [deg]

theta_vec = theta_mat(:,1);   % [180; 177; ...; 0] (61×1)
phi_vec   = phi_mat(1,:);     % [-180, -177, ..., 180] (1×121)

fprintf('데이터 로드 완료\n');
fprintf('  Theta: %.0f° ~ %.0f° (%.0f° 간격, %d개)\n', ...
    max(theta_vec), min(theta_vec), abs(theta_vec(2)-theta_vec(1)), nTheta);
fprintf('  Phi:   %.0f° ~ %.0f° (%.0f° 간격, %d개)\n', ...
    min(phi_vec), max(phi_vec), abs(phi_vec(2)-phi_vec(1)), nPhi);

%% ===== 2. 총 이득 계산 =====
% dBi → 선형 변환: G_lin = 10^(G_dBi / 10)
Gv_lin = 10.^(Gv_dBi_mat / 10);   % E_theta 선형 이득
Gh_lin = 10.^(Gh_dBi_mat / 10);   % E_phi   선형 이득

% 총 전력 이득 = 두 편광 성분의 합 (전력은 더해짐)
G_total_lin = Gv_lin + Gh_lin;
G_total_dBi = 10 * log10(G_total_lin);

[G_max_val, G_max_idx] = max(G_total_lin(:));
[r_max, c_max] = ind2sub(size(G_total_lin), G_max_idx);

fprintf('\n=== 이득 분석 ===\n');
fprintf('  최대 총 이득: %.4f dBi\n', G_total_dBi(r_max, c_max));
fprintf('  최대 위치: Theta = %.0f°, Phi = %.0f°\n', ...
    theta_vec(r_max), phi_vec(c_max));

%% ===== 3. 지향성(Directivity) 계산 =====
%
%  정의 (강의노트 4.1):
%    D(theta, phi) = 4*pi * U(theta, phi) / P_rad
%
%  여기서:
%    U(theta, phi) : 복사 강도 (radiation intensity) [W/sr]
%    P_rad = ∫∫ U(theta,phi) sin(theta) d_theta d_phi  [W]
%
%  측정 데이터에서 U(theta,phi) ∝ G_total_lin(theta,phi)
%  (far-field 에서 방사 패턴 형태가 곧 U의 형태)
%
%  수치 적분 (직사각형법):
%    P_rad ≈ Σ_i Σ_j  G_total_lin(i,j) * sin(theta_i) * dTheta * dPhi

dTheta_rad = abs(deg2rad(theta_vec(2) - theta_vec(1)));  % 3° → rad
dPhi_rad   = abs(deg2rad(phi_vec(2)   - phi_vec(1)));    % 3° → rad

sin_theta = sin(deg2rad(theta_vec));  % (61×1) sin 가중치

% 행렬 연산으로 적분: sin_theta를 열 방향으로 복제하여 element-wise 곱
P_rad = sum(sum(G_total_lin .* sin_theta)) * dTheta_rad * dPhi_rad;

% 최대 지향성
D_max_lin = 4 * pi * max(G_total_lin(:)) / P_rad;
D_max_dBi = 10 * log10(D_max_lin);

fprintf('\n=== 지향성(Directivity) 계산 결과 ===\n');
fprintf('  수치 적분 P_rad = %.6f\n', P_rad);
fprintf('  최대 U (선형)   = %.6f\n', max(G_total_lin(:)));
fprintf('  최대 지향성     = %.4f (선형) = %.2f dBi\n', D_max_lin, D_max_dBi);

%% ===== 4. 3D 방사 패턴 시각화 =====
%
%  구면좌표 → 직교좌표 변환:
%    x = r * sin(theta) * cos(phi)
%    y = r * sin(theta) * sin(phi)
%    z = r * cos(theta)
%
%  r = 정규화된 선형 이득 (0~1)

% theta, phi를 라디안으로 변환하고 meshgrid 생성
% phi_vec는 1×121 행벡터, theta_vec는 61×1 열벡터
[Phi_grid, Theta_grid] = meshgrid(deg2rad(phi_vec), deg2rad(theta_vec));
% 결과: Theta_grid, Phi_grid 모두 61×121

% 시각화용 정규화 (선형, 0~1)
R_norm = G_total_lin / max(G_total_lin(:));

X = R_norm .* sin(Theta_grid) .* cos(Phi_grid);
Y = R_norm .* sin(Theta_grid) .* sin(Phi_grid);
Z = R_norm .* cos(Theta_grid);

figure('Name', 'Project1 - 3D Radiation Pattern', 'Position', [50 50 850 650]);
surf(X, Y, Z, G_total_dBi, 'EdgeColor', 'none');
colorbar;
colormap jet;
xlabel('X'); ylabel('Y'); zlabel('Z');
title(sprintf('3D Radiation Pattern (Total Gain)  |  D_{max} = %.2f dBi', D_max_dBi), ...
    'FontSize', 13);
axis equal; grid on;
view(-45, 30);

%% ===== [선택] patternCustom 사용 시 (MATLAB Antenna Toolbox 필요) =====
%
%  patternCustom(magE, theta, phi)
%    - magE : [numel(theta) x numel(phi)] 행렬, 단위 dBi
%    - theta : 0~180 [deg], 증가 순서
%    - phi   : 0~360 또는 -180~180 [deg]
%
%  현재 데이터: theta가 180→0 순서이므로 flipud 필요

USE_PATTERN_CUSTOM = false;  % Antenna Toolbox 있으면 true로 변경

if USE_PATTERN_CUSTOM
    theta_asc = flipud(theta_vec);           % 0→180 (증가)
    G_flip    = flipud(G_total_dBi);         % theta 방향 뒤집기

    figure('Name', 'Project1 - patternCustom');
    patternCustom(G_flip, theta_asc, phi_vec);
    title(sprintf('3D Pattern (patternCustom)  |  D_{max} = %.2f dBi', D_max_dBi));
end

%% ===== 5. 2D 단면 패턴 (E-Plane & H-Plane) =====
%
%  E-Plane : 최대 복사 방향을 포함하는 전기장(E)의 면
%            주 편광(Vert=E_theta)의 최대 방향 → phi = 0°
%  H-Plane : E-Plane에 직교하는 자기장(H)의 면 → phi = 90°

phi0_idx  = find(phi_vec == 0);
phi90_idx = find(phi_vec == 90);

figure('Name', 'Project1 - 2D Cut Patterns', 'Position', [50 50 1100 480]);

% --- E-Plane (phi = 0°) ---
subplot(1, 2, 1);
plot(theta_vec, G_total_dBi(:, phi0_idx), 'b-',  'LineWidth', 2.5); hold on;
plot(theta_vec, Gv_dBi_mat(:, phi0_idx),  'r--', 'LineWidth', 1.5);
plot(theta_vec, Gh_dBi_mat(:, phi0_idx),  'g--', 'LineWidth', 1.5);
grid on;
xlabel('Theta (deg)', 'FontSize', 11);
ylabel('Gain (dBi)', 'FontSize', 11);
title('E-Plane (φ = 0°)', 'FontSize', 12);
legend('Total', 'Vert (E_\theta)', 'Horiz (E_\phi)', 'Location', 'best');
xlim([0 180]);
set(gca, 'XDir', 'reverse');  % theta=0이 오른쪽 (z축 방향)

% --- H-Plane (phi = 90°) ---
subplot(1, 2, 2);
plot(theta_vec, G_total_dBi(:, phi90_idx), 'b-',  'LineWidth', 2.5); hold on;
plot(theta_vec, Gv_dBi_mat(:, phi90_idx),  'r--', 'LineWidth', 1.5);
plot(theta_vec, Gh_dBi_mat(:, phi90_idx),  'g--', 'LineWidth', 1.5);
grid on;
xlabel('Theta (deg)', 'FontSize', 11);
ylabel('Gain (dBi)', 'FontSize', 11);
title('H-Plane (φ = 90°)', 'FontSize', 12);
legend('Total', 'Vert (E_\theta)', 'Horiz (E_\phi)', 'Location', 'best');
xlim([0 180]);
set(gca, 'XDir', 'reverse');

sgtitle(sprintf('2D Pattern Cuts  |  D_{max} = %.2f dBi', D_max_dBi), 'FontSize', 13);

%% ===== 6. 2D 패턴: -180°~180° 표현 (E-Plane & H-Plane) =====
%
%  구성 원리:
%    angle < 0 : phi=180°(E) 또는 phi=-90°(H) 쪽 → theta 180→0
%    angle ≥ 0 : phi=0°(E)  또는 phi=90°(H)  쪽 → theta 0→180
%
%    E-Plane: xz 평면 (phi=0°/180°)
%    H-Plane: yz 평면 (phi=90°/-90°)

phi180_idx = find(phi_vec == 180);
phim90_idx = find(phi_vec == -90);

angle_axis = -180:3:180;   % 121점

% ── E-Plane ──
% 왼쪽(-180→0): phi=180°, theta 내림차순(180→0) → 그대로 사용
G_E_left  = G_total_dBi(:, phi180_idx);                      % 61점
Gv_E_left = Gv_dBi_mat(:, phi180_idx);
Gh_E_left = Gh_dBi_mat(:, phi180_idx);
% 오른쪽(3→180): phi=0°, theta 오름차순(3→180) → flipud 후 끝점 제외
G_E_right  = flipud(G_total_dBi(1:60, phi0_idx));            % 60점
Gv_E_right = flipud(Gv_dBi_mat(1:60, phi0_idx));
Gh_E_right = flipud(Gh_dBi_mat(1:60, phi0_idx));

G_E_full  = [G_E_left;  G_E_right];   % 121점
Gv_E_full = [Gv_E_left; Gv_E_right];
Gh_E_full = [Gh_E_left; Gh_E_right];

% ── H-Plane ──
G_H_left  = G_total_dBi(:, phim90_idx);
Gv_H_left = Gv_dBi_mat(:, phim90_idx);
Gh_H_left = Gh_dBi_mat(:, phim90_idx);

G_H_right  = flipud(G_total_dBi(1:60, phi90_idx));
Gv_H_right = flipud(Gv_dBi_mat(1:60, phi90_idx));
Gh_H_right = flipud(Gh_dBi_mat(1:60, phi90_idx));

G_H_full  = [G_H_left;  G_H_right];
Gv_H_full = [Gv_H_left; Gv_H_right];
Gh_H_full = [Gh_H_left; Gh_H_right];

% ── 플롯 ──
figure('Name','Project1 - 2D Patterns (-180~180)', 'Position',[50 50 1100 480]);

subplot(1,2,1);
plot(angle_axis, G_E_full,  'b-',  'LineWidth', 2.5); hold on;
plot(angle_axis, Gv_E_full, 'r--', 'LineWidth', 1.5);
plot(angle_axis, Gh_E_full, 'g--', 'LineWidth', 1.5);
grid on;
xlabel('Angle (deg)', 'FontSize', 11);
ylabel('Gain (dBi)', 'FontSize', 11);
title('E-Plane  (φ=0°/180°)', 'FontSize', 12);
legend('Total', 'Vert (E_\theta)', 'Horiz (E_\phi)', 'Location', 'best');
xlim([-180 180]);
xticks(-180:30:180);
xline(0, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');

subplot(1,2,2);
plot(angle_axis, G_H_full,  'b-',  'LineWidth', 2.5); hold on;
plot(angle_axis, Gv_H_full, 'r--', 'LineWidth', 1.5);
plot(angle_axis, Gh_H_full, 'g--', 'LineWidth', 1.5);
grid on;
xlabel('Angle (deg)', 'FontSize', 11);
ylabel('Gain (dBi)', 'FontSize', 11);
title('H-Plane  (φ=90°/-90°)', 'FontSize', 12);
legend('Total', 'Vert (E_\theta)', 'Horiz (E_\phi)', 'Location', 'best');
xlim([-180 180]);
xticks(-180:30:180);
xline(0, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');

sgtitle(sprintf('2D Pattern Cuts (-180°~180°)  |  D_{max} = %.2f dBi', D_max_dBi), ...
    'FontSize', 13);

%% ===== 7. 결과 요약 출력 =====
fprintf('\n========== Project 1 결과 요약 ==========\n');
fprintf('  데이터 크기    : %d × %d (theta × phi)\n', nTheta, nPhi);
fprintf('  Theta 범위     : 0° ~ 180° (3° 간격)\n');
fprintf('  Phi 범위       : -180° ~ 180° (3° 간격)\n');
fprintf('  최대 총 이득   : %.4f dBi  @ Theta=%.0f°, Phi=%.0f°\n', ...
    G_total_dBi(r_max, c_max), theta_vec(r_max), phi_vec(c_max));
fprintf('  최대 지향성    : %.4f (선형) = %.2f dBi\n', D_max_lin, D_max_dBi);
fprintf('==========================================\n');
