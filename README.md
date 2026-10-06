# OAI-Based Interference Suppression

OpenAirInterface(OAI) 5G 기지국에 **소프트웨어 폴딩, 신호 복원 및 두 안테나 기반 공간 필터링**을 통합한 연구용 프로그램입니다. USRP에서 수신한 상향링크 I/Q 신호를 OAI의 FFT 이전 단계에서 처리하여, 기존 기지국 수신 체계와 연동합니다.

기반 소스는 OAI `v2.3.0` 체크아웃 시점의 커밋 `8bf6d5d7da`입니다. 아래 설명은 이 저장소에서 추가한 기능을 대상으로 하며, OAI 원본 소개와 문서는 하단에 유지합니다.

## 수신 처리 구조

```text
USRP 수신 I/Q 데이터 (ANT0, ANT1)
  → 상향링크·혼합 슬롯의 OFDM 심볼별 처리
  → I/Q 정규화 및 모듈로 기반 소프트웨어 폴딩
  → 차분·누적합 기반 신호 복원
  → 원본 첫 샘플을 기준으로 복원 오프셋 보정
  → 두 안테나의 상관관계에 기반한 공간 필터링
  → ANT0에 결과 저장, ANT1은 0으로 설정
  → OAI FFT 및 후속 물리계층 수신 처리
```

심볼별 처리는 OpenMP로 병렬화하며, 초기화 시 할당한 스레드별 작업 버퍼를 재사용합니다. 공간 필터링은 두 수신 신호의 에너지와 복소 상관값으로 정규화된 결합 가중치를 계산하여 간섭 억제를 수행합니다.

## 자체 개발 및 변경 부분

| 파일 | 역할 |
|---|---|
| [reconstruction_test_multi.c](executables/reconstruction_test_multi.c) | 현재 호출되는 폴딩·복원·오프셋 보정 및 공간 필터링 구현 |
| [reconstruction_test_multi.h](executables/reconstruction_test_multi.h) | 스레드별 작업 버퍼 및 처리 크기 정의 |
| [reconstruction_test.c](executables/reconstruction_test.c) | 초기 복원 실험 구현. 빌드에 포함되지만 현재 수신 경로에서 호출하지 않음 |
| [nr-ru.c](executables/nr-ru.c) | FFT 이전 처리 삽입, OFDM 심볼 병렬 처리, 수신 버퍼 추적 |
| [CMakeLists.txt](CMakeLists.txt) | 복원 모듈 빌드 및 OpenMP 연결 |
| [T_messages.txt](common/utils/T/T_messages.txt) / [usrp_lib.cpp](radio/USRP/usrp_lib.cpp) | 안테나별 추적 이벤트 및 데이터 수집 위치 변경 |
| [USRP N310 실험 설정](targets/PROJECTS/GENERIC-NR-5GC/CONF/gnb.band78.sa.fr1.106PRB.usrpn310_mod.conf) | Band 78, 106 PRB, 송신 1채널·수신 2채널 실험 구성 |

## 실행 환경 및 빌드

Linux PC, UHD 지원 USRP, 두 수신 채널, OAI 빌드 의존성과 OpenMP를 지원하는 C/C++ 컴파일러가 필요합니다. 단말 접속 및 상향링크 전송 시험에는 별도의 5G 코어망과 시험용 단말도 필요합니다. 의존성 설치 절차는 [OAI 빌드 문서](doc/BUILD.md)를 참고하세요.

```bash
git clone https://github.com/chkim04/oai_based_interference_suppression.git
cd oai_based_interference_suppression/cmake_targets
./build_oai -w USRP --gNB
```

빌드 로그에서 OpenMP가 검출되고 활성화되는지 확인합니다. 현재 추가 모듈은 OpenMP 헤더와 런타임 함수를 직접 사용합니다.

## 설정 및 실행 방법

1. USRP와 PC의 네트워크 연결 및 두 수신 채널을 준비합니다.
2. 위 실험 설정 파일의 `sdr_addrs`, 주파수, 대역폭, 수신 이득, PLMN 및 코어망 연결 정보를 실제 환경에 맞게 수정합니다. 파일에는 장치 주소 `192.168.10.2` / `192.168.11.2`와 외부 클록·시간 동기원이 지정되어 있으므로 장치 구성에 맞춰 변경해야 합니다.
3. 복원 조건을 변경하려면 `reconstruction_test_multi.c`의 `LAMBDA`(현재 `0.05`)와 `RECON_N`(현재 `2`)을 수정한 후 다시 빌드합니다. 현재 이 값들은 명령행 옵션이 아닙니다.
4. 빌드 결과 디렉터리에서 기지국을 실행합니다. 아래 스레드 수 `4`는 실행 예시이며, CPU와 처리 시간에 맞게 조정합니다. 현재 `MAX_THREADS`는 `32`, `MAX_SYM_LEN`은 CP를 포함하여 `4096` 샘플이므로 이 범위 안에서 운용합니다.

```bash
# 저장소 루트에서 실행
cd cmake_targets/ran_build/build
sudo env OMP_NUM_THREADS=4 ./nr-softmodem \
  -O ../../../targets/PROJECTS/GENERIC-NR-5GC/CONF/gnb.band78.sa.fr1.106PRB.usrpn310_mod.conf
```

5. USRP 초기화와 코어망 연결을 확인하고 단말을 접속시킨 뒤 상향링크 데이터를 전송합니다. 수신 경로에 통합된 복원 및 공간 필터링이 자동으로 실행됩니다.
6. 실행 로그에서 단말 접속과 상향링크 수신 상태를 확인합니다. 시험 완료 후 단말 트래픽을 중지하고 기지국을 종료합니다. 위 명령은 기본 실행 형식이며, 환경별 추가 옵션은 실제 운용 설정에 맞춥니다.

## 데이터 수집

OAI T tracer의 `T_USRP_RX_ANT0`, `T_USRP_RX_ANT1` 이벤트로 RU 수신 버퍼를 추적할 수 있습니다. 사용법은 [T tracer 문서](common/utils/T/DOC/T.md)와 [기록 방법](common/utils/T/DOC/T/record.md)을 참고하세요.

현재 추적 위치는 RU의 해당 슬롯 처리 이후입니다. 공간 필터링을 거친 구간의 ANT0에는 처리 결과가, ANT1에는 0이 기록되므로 두 안테나의 원시 입력을 보존한 덤프로 해석하면 안 됩니다. 로컬 실험 캡처 파일 `common/utils/T/tracer/tp_uplink_rx_data_wireless_*`는 Git 추적에서 제외됩니다.

## 현재 구현 범위

- 폴딩은 USRP가 디지털화한 신호에 소프트웨어로 적용합니다. 실제 폴딩 ADC 하드웨어 입력을 직접 받는 구현은 아닙니다.
- 복원 오프셋 보정은 폴딩 이전 입력의 첫 샘플을 참조합니다. 폴딩된 데이터만으로 수행하는 완전한 블라인드 복원으로 해석하지 않습니다.
- 주파수 회전·역회전 함수는 구현되어 있지만 현재 수신 경로의 호출은 비활성화되어 있습니다. 주석 처리된 노치 필터도 현재 동작에 포함되지 않습니다.
- 현재 실험 경로는 두 수신 채널을 기준으로 구성되어 있습니다. 단일 수신 채널을 사용하려면 추적 코드 등의 채널 접근도 함께 검토해야 합니다.
- 간섭 억제 성능과 처리 시간은 입력 신호, 채널 및 장비 설정에 따라 검증해야 합니다. 이 README에는 새로 측정한 성능 수치나 검증되지 않은 성능 보장을 포함하지 않습니다.

## OAI 원본 및 라이선스

기존 OAI 소스의 출처와 라이선스 정보는 [LICENSE](LICENSE), [NOTICE.md](NOTICE.md)에 있습니다. 다음은 OAI 원본 프로젝트의 안내입니다.

---

<h1 align="center">
    <a href="https://openairinterface.org/"><img src="https://openairinterface.org/wp-content/uploads/2015/06/cropped-oai_final_logo.png" alt="OAI" width="550"></a>
</h1>

<p align="center">
    <a href="https://gitlab.eurecom.fr/oai/openairinterface5g/-/blob/master/LICENSE"><img src="https://img.shields.io/badge/license-OAI--Public--V1.1-blue" alt="License"></a>
    <a href="https://releases.ubuntu.com/22.04/"><img src="https://img.shields.io/badge/OS-Ubuntu22-Green" alt="Supported OS Ubuntu 22"></a>
    <a href="https://releases.ubuntu.com/24.04/"><img src="https://img.shields.io/badge/OS-Ubuntu24-Green" alt="Supported OS Ubuntu 24"></a>
    <a href="https://www.redhat.com/en/technologies/linux-platforms/enterprise-linux"><img src="https://img.shields.io/badge/OS-RHEL9-Green" alt="Supported OS RHEL9"></a>
    <a href="https://getfedora.org/en/workstation/"><img src="https://img.shields.io/badge/OS-Fedore41-Green" alt="Supported OS Fedora 41"></a>
</p>

<p align="center">
    <a href="https://gitlab.eurecom.fr/oai/openairinterface5g/-/releases"><img alt="GitLab Release (custom instance)" src="https://img.shields.io/gitlab/v/release/oai/openairinterface5g?gitlab_url=https%3A%2F%2Fgitlab.eurecom.fr&include_prereleases&sort=semver"></a>
</p>

<p align="center">
    <a href="https://jenkins-oai.eurecom.fr/job/RAN-Ubuntu18-Image-Builder/"><img src="https://img.shields.io/jenkins/build?jobUrl=https%3A%2F%2Fjenkins-oai.eurecom.fr%2Fjob%2FRAN-Ubuntu18-Image-Builder%2F&label=build-Ubuntu-x86%20Images"></a>
    <a href="https://jenkins-oai.eurecom.fr/job/RAN-RHEL8-Cluster-Image-Builder/"><img src="https://img.shields.io/jenkins/build?jobUrl=https%3A%2F%2Fjenkins-oai.eurecom.fr%2Fjob%2FRAN-RHEL8-Cluster-Image-Builder%2F&label=build-UBI-x86%20Images"></a>
    <a href="https://jenkins-oai.eurecom.fr/job/RAN-Ubuntu-ARM-Image-Builder/"><img src="https://img.shields.io/jenkins/build?jobUrl=https%3A%2F%2Fjenkins-oai.eurecom.fr%2Fjob%2FRAN-Ubuntu-ARM-Image-Builder%2F&label=build-Ubuntu-ARM%20Images"></a>
</p>

<p align="center">
  <a href="https://hub.docker.com/r/oaisoftwarealliance/oai-gnb"><img alt="Docker Pulls" src="https://img.shields.io/docker/pulls/oaisoftwarealliance/oai-gnb?label=gNB%20docker%20pulls"></a>
  <a href="https://hub.docker.com/r/oaisoftwarealliance/oai-nr-ue"><img alt="Docker Pulls" src="https://img.shields.io/docker/pulls/oaisoftwarealliance/oai-nr-ue?label=NR-UE%20docker%20pulls"></a>
  <a href="https://hub.docker.com/r/oaisoftwarealliance/oai-enb"><img alt="Docker Pulls" src="https://img.shields.io/docker/pulls/oaisoftwarealliance/oai-enb?label=eNB%20docker%20pulls"></a>
  <a href="https://hub.docker.com/r/oaisoftwarealliance/oai-lte-ue"><img alt="Docker Pulls" src="https://img.shields.io/docker/pulls/oaisoftwarealliance/oai-lte-ue?label=LTE-UE%20docker%20pulls"></a>
  <a href="https://hub.docker.com/r/oaisoftwarealliance/oai-nr-cuup"><img alt="Docker Pulls" src="https://img.shields.io/docker/pulls/oaisoftwarealliance/oai-nr-cuup?label=NR-CUUP%20docker%20pulls"></a>
</p>

# OpenAirInterface License #

 *  [OAI License Model](http://www.openairinterface.org/?page_id=101)
 *  [OAI License v1.1 on our website](http://www.openairinterface.org/?page_id=698)

It is distributed under **OAI Public License V1.1**.

The license information is distributed under [LICENSE](LICENSE) file in the same directory.

Please see [NOTICE](NOTICE.md) file for third party software that is included in the sources.

# Where to Start #

 *  [General overview of documentation](./doc/README.md)
 *  [The implemented features](./doc/FEATURE_SET.md)
 *  [System Requirements for Using OAI Stack](./doc/system_requirements.md)
 *  [How to build](./doc/BUILD.md)
 *  [How to run the modems](./doc/RUNMODEM.md)

Not all information is available in a central place, and information for
specific sub-systems might be available in the corresponding sub-directories.
To find all READMEs, this command might be handy:

```
find . -iname "readme*"
```

# RAN repository structure #

The OpenAirInterface (OAI) software is composed of the following parts: 

```
openairinterface5g
├── charts
├── ci-scripts        : Meta-scripts used by the OSA CI process. Contains also configuration files used day-to-day by CI.
├── CMakeLists.txt    : Top-level CMakeLists.txt for building
├── cmake_targets     : Build utilities to compile (simulation, emulation and real-time platforms), and generated build files.
├── common            : Some common OAI utilities, some other tools can be found at openair2/UTILS.
├── doc               : Documentation
├── docker            : Dockerfiles to build for Ubuntu and RHEL
├── executables       : Top-level executable source files (gNB, eNB, ...)
├── maketags          : Script to generate emacs tags.
├── nfapi             : (n)FAPI code for MAC-PHY interface
├── openair1          : Layer 1 (3GPP LTE Rel-10/12 PHY, NR Rel-15 PHY)
├── openair2          : Layer 2 (3GPP LTE Rel-10 MAC/RLC/PDCP/RRC/X2AP, LTE Rel-14 M2AP, NR Rel-15+ MAC/RLC/PDCP/SDAP/RRC/X2AP/F1AP/E1AP), E2AP
├── openair3          : Layer 3 (3GPP LTE Rel-10 S1AP/GTP, NR Rel-15 NGAP/GTP)
├── openshift         : OpenShift helm charts for some deployment options of OAI
├── radio             : Drivers for various radios such as USRP, AW2S, RFsim, 7.2 FHI, ...
├── targets           : Some configuration files; only historical relevance, and might be deleted in the future
└── tools             : Tools for use by the developers/ci machines: code analysis and formatting
```

# How to get support from the OAI Community # 

You can ask your question on the [mailing lists](https://gitlab.eurecom.fr/oai/openairinterface5g/-/wikis/MailingList).

Your email should contain below information:

- A clear subject in your email.
- For all the queries there should be [Query\] in the subject of the email and for problems there should be [Problem\].
- In case of a problem, add a small description.
- Do not share any photos unless you want to share a diagram.
- OAI gNB/DU/CU/CU-CP/CU-UP configuration file in `.conf` format only.
- Logs of OAI gNB/DU/CU/CU-CP/CU-UP in `.log` or `.txt` format only.
- In case your question is related to performance, include a small description of the machine (Operating System, Kernel version, CPU, RAM and networking card) and diagram of your testing environment.
- Known/open issues are present on [GitLab](https://gitlab.eurecom.fr/oai/openairinterface5g/-/issues), so keep checking.

Always remember a structured email will help us understand your issues quickly.
