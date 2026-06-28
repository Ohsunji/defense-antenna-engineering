%% Project 2: CP 안테나 측정데이터를 사용한 성능 추출
%  이론: 강의노트 6 (CP Antenna Measurement)
%
%  입력 데이터 (A~F 열):
%    A: Theta [deg]
%    B: Phi   [deg]
%    C: Vert  Gain  (E_V 성분, 수직 편광) [dBi]
%    D: Vert  Phase (E_V 성분) [deg]
%    E: Horiz Gain  (E_H 성분, 수평 편광) [dBi]
%    F: Horiz Phase (E_H 성분) [deg]
%
%  계산 목표 (G~I 열):
%    G: LHCP Gain [dBi]
%    H: RHCP Gain [dBi]
%    I: Axial Ratio [dB]

clear; clc;

%% ===== 데이터 입력 =====
% 형식: [Theta, Phi, Gv(dBi), Pv(deg), Gh(dBi), Ph(deg)]
% 아래에 실제 측정 데이터 입력 (슬라이드 표의 값)
data = [
     0,  0,  -4.89, 359.48,  -9.48, 258.05;
     2,  0,  -4.89, 359.64,  -9.44, 258.40;
     4,  0,  -4.89, 359.54,  -9.40, 258.47;
     6,  0,  -4.89, 359.42,  -9.37, 258.42;
     8,  0,  -4.90, 359.35,  -9.36, 258.36;
    10,  0,  -4.91, 359.33,  -9.34, 258.21;
    12,  0,  -4.93, 359.26,  -9.32, 257.98;
    14,  0,  -4.96, 359.11,  -9.31, 257.60;
    16,  0,  -4.98, 359.02,  -9.29, 257.26;
    18,  0,  -5.02, 359.03,  -9.28, 257.09;
    20,  0,  -5.06, 358.86,  -9.29, 256.82;
    % 데이터가 더 있으면 여기에 추가
];

Theta   = data(:,1);
Phi     = data(:,2);
Gv_dBi  = data(:,3);   % Vert Gain [dBi]
Pv_deg  = data(:,4);   % Vert Phase [deg]
Gh_dBi  = data(:,5);   % Horiz Gain [dBi]
Ph_deg  = data(:,6);   % Horiz Phase [deg]

N = size(data, 1);
fprintf('데이터 로드 완료: %d개 포인트\n', N);

%% ===== 이론 설명 =====
%
%  CP 분해 공식 (강의노트 6, Ludwig 3rd definition 기반):
%
%  수직/수평 선형 편광 성분을 복소 페이저(complex phasor)로 표현:
%    E_V = |E_V| * exp(j * phi_V)    |E_V| = 10^(Gv/20)  [field amplitude]
%    E_H = |E_H| * exp(j * phi_H)    |E_H| = 10^(Gh/20)
%
%  LHCP / RHCP 성분 분해 (강의노트 6, p.45 수식 기준):
%    rho_RHCP = (v_hat + j*h_hat) / sqrt(2)  →  E_RHCP = (E_V + j * E_H) / sqrt(2)
%    rho_LHCP = (v_hat - j*h_hat) / sqrt(2)  →  E_LHCP = (E_V - j * E_H) / sqrt(2)
%
%  |E_RHCP|^2 전개:
%    = (|E_V|^2 + |E_H|^2 - 2|E_V||E_H| * sin(phi_H - phi_V)) / 2
%
%  |E_LHCP|^2 전개:
%    = (|E_V|^2 + |E_H|^2 + 2|E_V||E_H| * sin(phi_H - phi_V)) / 2
%
%  이득 변환:
%    G_LHCP [dBi] = 10 * log10(|E_LHCP|^2)
%    G_RHCP [dBi] = 10 * log10(|E_RHCP|^2)
%
%  축비 (Axial Ratio):
%    AR_linear = (|E_LHCP| + |E_RHCP|) / | |E_LHCP| - |E_RHCP| |
%    AR [dB]   = 20 * log10(AR_linear)
%    (완벽한 CP: AR = 0 dB / 선형 편광: AR → ∞ dB)

%% ===== Step 1: 선형 이득 및 위상 변환 =====
% 전계 진폭 (field amplitude): |E| = 10^(G_dBi / 20)
%   (전력 이득 G_lin = |E|^2 = 10^(G_dBi/10) 이므로 |E| = 10^(G_dBi/20))
Ev_amp = 10.^(Gv_dBi / 20);   % |E_V| (field amplitude)
Eh_amp = 10.^(Gh_dBi / 20);   % |E_H| (field amplitude)

% 위상을 라디안으로 변환
Pv_rad = deg2rad(Pv_deg);
Ph_rad = deg2rad(Ph_deg);

%% ===== Step 2: 복소 페이저 구성 =====
E_V = Ev_amp .* exp(1j * Pv_rad);   % E_V 복소 페이저
E_H = Eh_amp .* exp(1j * Ph_rad);   % E_H 복소 페이저

%% ===== Step 3: LHCP / RHCP 분해 =====
E_RHCP = (E_V + 1j * E_H) / sqrt(2);   % rho_RHCP = (v + jh)/sqrt(2)
E_LHCP = (E_V - 1j * E_H) / sqrt(2);   % rho_LHCP = (v - jh)/sqrt(2)

% 전력 (|E|^2)
P_LHCP = abs(E_LHCP).^2;
P_RHCP = abs(E_RHCP).^2;

% dBi 변환
G_LHCP_dBi = 10 * log10(P_LHCP);   % 열 G
G_RHCP_dBi = 10 * log10(P_RHCP);   % 열 H

%% ===== Step 4: 축비(Axial Ratio) 계산 =====
% 전계 진폭 (CP 성분)
E_L_amp = abs(E_LHCP);   % |E_LHCP|
E_R_amp = abs(E_RHCP);   % |E_RHCP|

% AR = (|E_L| + |E_R|) / | |E_L| - |E_R| |
AR_linear = (E_L_amp + E_R_amp) ./ abs(E_L_amp - E_R_amp);
AR_dB = 20 * log10(AR_linear);    % 열 I

%% ===== 결과 출력 =====
fprintf('\n========== Project 2 계산 결과 ==========\n');
fprintf('%-6s %-5s | %-8s %-8s | %-10s %-10s | %-10s\n', ...
    'Theta', 'Phi', 'Gv(dBi)', 'Gh(dBi)', 'LHCP(dBi)', 'RHCP(dBi)', 'AR(dB)');
fprintf('%s\n', repmat('-', 1, 75));

for i = 1:N
    fprintf('%-6.0f %-5.0f | %-8.2f %-8.2f | %-10.2f %-10.2f | %-10.2f\n', ...
        Theta(i), Phi(i), Gv_dBi(i), Gh_dBi(i), ...
        G_LHCP_dBi(i), G_RHCP_dBi(i), AR_dB(i));
end

%% ===== 검증: 전력 보존 확인 =====
% LHCP + RHCP 전력 = 전체 전력 이어야 함
fprintf('\n===== 검증: 전력 보존 (LHCP + RHCP = Vtotal + Htotal) =====\n');
P_total_input = 10.^(Gv_dBi/10) + 10.^(Gh_dBi/10);
P_total_output = P_LHCP + P_RHCP;
fprintf('최대 오차: %.2e (거의 0이면 정상)\n', max(abs(P_total_input - P_total_output)));

%% ===== 시각화 =====
figure('Name', 'Project2 - CP Performance', 'Position', [50 50 1100 500]);

% --- subplot 1: LHCP / RHCP 이득 ---
subplot(1, 2, 1);
plot(Theta, G_LHCP_dBi, 'b-o', 'LineWidth', 2, 'MarkerSize', 5); hold on;
plot(Theta, G_RHCP_dBi, 'r-s', 'LineWidth', 2, 'MarkerSize', 5);
plot(Theta, Gv_dBi, 'k--', 'LineWidth', 1.2);
plot(Theta, Gh_dBi, 'k:', 'LineWidth', 1.2);
grid on;
xlabel('Theta (deg)', 'FontSize', 11);
ylabel('Gain (dBi)', 'FontSize', 11);
title('LHCP / RHCP Gain', 'FontSize', 12);
legend('LHCP (G열)', 'RHCP (H열)', 'Vert (입력)', 'Horiz (입력)', 'Location', 'best');
xlim([min(Theta) max(Theta)]);

% --- subplot 2: 축비 ---
subplot(1, 2, 2);
plot(Theta, AR_dB, 'm-^', 'LineWidth', 2, 'MarkerSize', 6);
hold on;
yline(3, 'k--', '3 dB 기준선', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
grid on;
xlabel('Theta (deg)', 'FontSize', 11);
ylabel('Axial Ratio (dB)', 'FontSize', 11);
title('축비 (Axial Ratio)', 'FontSize', 12);
legend('AR (I열)', 'Location', 'best');
xlim([min(Theta) max(Theta)]);

sgtitle('Project 2: CP 안테나 성능 추출 결과', 'FontSize', 13);

%% ===== Excel 출력 (선택) =====
% 결과를 Excel로 저장하려면 아래 주석 해제
% result_table = table(Theta, Phi, Gv_dBi, Pv_deg, Gh_dBi, Ph_deg, ...
%     G_LHCP_dBi, G_RHCP_dBi, AR_dB, ...
%     'VariableNames', {'Theta','Phi','Gv_dBi','Pv_deg','Gh_dBi','Ph_deg', ...
%                       'LHCP_dBi','RHCP_dBi','AR_dB'});
% writetable(result_table, 'Project2_result.xlsx');
% fprintf('Excel 저장 완료: Project2_result.xlsx\n');
