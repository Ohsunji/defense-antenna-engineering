# 국방안테나공학 프로젝트 결과 정리

## Project 1 - 안테나 원거리장 패턴 및 링크버짓

### 내용
- 반파장 다이폴 안테나의 이론적 방사 패턴 계산 및 3D 시각화
- ULA(균일 선형 배열) 빔 패턴 분석 (배열 요소 수, 간격 변화)
- Friis 전송 방정식 기반 링크버짓: 주파수·거리·안테나 이득 파라미터 스윕

### 주요 결과
| 항목 | 결과 |
|------|------|
| 다이폴 최대 이득 | ~2.15 dBi |
| 4소자 ULA 메인빔 폭 | ~26° (d=λ/2) |
| 링크버짓 마진 (10 km, 1 GHz) | +12.3 dB |

---

## Project 2 - 원형편파 성능 분석 및 근접장 측정

### 내용
- 평면 근접장 스캔 데이터(SGH)로부터 PWS(평면파 스펙트럼) 변환
- Back-Projection 알고리즘으로 원거리장 재구성
- 원형편파(CP) 성능: AR(Axial Ratio), 편파 효율 계산

### 주요 결과
| 항목 | 결과 |
|------|------|
| 재구성 원거리장 vs 직접 측정 오차 | < 1 dB (주빔 영역) |
| CP Axial Ratio @ 주빔 | ~1.2 dB |
| 편파 효율 | ~98% |

---

## Project 3 - 빔포밍 배열안테나

### 내용
- ULA 빔포밍 가중치 설계: 균일(Uniform), Chebyshev, Taylor 윈도우
- 부엽 레벨(SLL) vs 메인빔 폭 트레이드오프 분석
- 멀티빔 생성: 위상 천이로 스캔 각도 제어

### 주요 결과
| 가중치 | SLL | HPBW |
|--------|-----|------|
| Uniform | -13 dB | 좁음 |
| Chebyshev (-30dB) | -30 dB | 넓어짐 |
| Taylor (n̄=5) | -25 dB | 중간 |

---

## 기말 - Infinitesimal Dipole Modeling (IDM)

### 내용
- IDM 기반 근접장 측정 보정 알고리즘 구현
- 다이폴 모델링 파라미터(위치, 방향, 크기)를 최적화하여 근접장 재현
- Back-Projection 결과와 IDM 결과 비교

### 주요 결과
- IDM 재구성 원거리장 패턴 vs 측정값: RMSE < 0.8 dB
- 계산 시간: 해석적 IDM이 Full-wave 대비 ~200× 빠름

---

## 파일 구조

```
src/
├── project1_antenna_pattern.m      # 방사 패턴 + 3D 시각화
├── project1_link_budget.m          # Friis 링크버짓 계산
├── project2_far_field.m            # PWS 변환 + 원거리장 재구성
├── project2_CP_performance.m       # 원형편파 성능 분석
├── project3_beamforming.m          # 빔포밍 가중치 설계
├── advanced_dipole_modeling.m      # IDM 기반 근접장 보정
└── BackProjection_PlanarNearField.m # Back-Projection 알고리즘
```

## 개발환경
- MATLAB R2024b
- Antenna Toolbox
