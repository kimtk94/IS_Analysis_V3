# MasterOmics 공통 데이터 수집기

공개 파일 다운로드, 중단 후 재개, SHA256 검증, 승인된 로컬 파일 staging을 같은 코어로 처리한다. 실제 생산 DB를 이번 구현에서 다운로드하지 않았다. DB별 column 변환과 분석 실행은 다음 단계다.

## 실행 조건

`resources plan`은 오프라인 명세 검사이고 `resources collect`만 실제 파일을 수집한다. 반드시 명시적 ID 선택, exact artifact, release, SHA256, license 검토가 필요하다. 유전 자료는 genome build와 ancestry도 고정한다. 포털 URL 대신 정확한 파일 URL을 지정한다. URL은 credentials가 없는 HTTPS만 허용한다. 토큰·쿠키·가입·승인 절차는 구현하지 않는다.

공개 remote 파일의 binding은 다음과 같다. 아래 null은 placeholder이므로 실행 시 보류된다. checksum은 실제 확인된 파일 또는 신뢰할 수 있는 배포 출처의 digest로 고정한다.

```json
{
  "ukb_ppp_pqtl": {
    "artifact_url": null,
    "file_name": null,
    "release": null,
    "sha256": null,
    "expected_bytes": null,
    "genome_build": null,
    "ancestry": null,
    "license_reviewed": false
  }
}
```

`expected_bytes`는 선택적 양의 정수다. 기록할 때 압축 파일이면 **압축 상태의 실제 바이트 수**를 사용한다. SHA256도 동일한 파일 바이트에 대한 값이다. 인위적인 예시 digest로 생산 데이터를 받지 않는다.

서버에서 기존 branch를 업데이트한 뒤 저장소 루트에서 실행한다. 예시 bindings 경로는 사용자가 검토한 로컬 명세로 교체한다.

```bash
python3 -S -m masteromics resources plan \
  --ids ukb_ppp_pqtl \
  --bindings /srv/is-analysis/config/masteromics/resources.bindings.json

bash server/masteromics_collect.sh \
  --ids ukb_ppp_pqtl \
  --bindings /srv/is-analysis/config/masteromics/resources.bindings.json \
  --root /srv/is-analysis/data/masteromics \
  --retries 2 --timeout 60
```

`MASTEROMICS_PYTHON`을 지정할 수도 있으나 수집기는 Python 표준 라이브러리만 필요하다. 설치를 실행하지 않는다. 서버 wrapper는 `set -e` 없이 종료 코드를 전달한다.

## 승인된 로컬 파일

등록/통제 접근 소스는 remote URL로 수집하지 않는다. 사용자가 이미 적법하게 확보한 파일만 `local_path`로 staging한다. `access_authorized=true`와 비밀정보를 포함하지 않는 `authorization_ref`가 함께 있어야 한다. 권한 승인 여부 자체를 자동 검증하거나 새로 발급하는 기능은 없다. reference는 내부 승인/사용계약 기록을 찾기 위한 식별자이며 토큰·비밀번호·승인 문서 본문을 넣지 않는다.

```json
{
  "koges_cohort": {
    "local_path": "/absolute/path/to/approved/file.tsv",
    "file_name": "cohort.tsv",
    "release": null,
    "sha256": null,
    "genome_build": null,
    "ancestry": null,
    "license_reviewed": false,
    "access_authorized": false,
    "authorization_ref": null
  }
}
```

local_path는 절대 경로여야 한다. `artifact_url`은 로컬 모드에서 사용하지 않으며 소스의 `access_type=controlled` 표시는 그대로 남는다. 공개 로컬 파일은 license/checksum 등 기본 조건을 충족하면 권한 reference 없이 staging할 수 있다. 원본을 변경하거나 이동하지 않고 별도의 수집 root로 복사한다.

## 저장 및 재개

| 출력 | 의미 |
|---|---|
| `raw/<resource_id>/<file_name>` | 검증 후 승격된 원본 바이트 파일 |
| `<file_name>.partial`, `<file_name>.partial.json` | 재개용 중간 바이트와 명세 identity/ETag |
| `<file_name>.acquisition.json` | 최신 시도의 status, identity, config/code SHA256, 시작·종료 시각, 오류 |
| `<file_name>.history.jsonl` | 잠금 아래 추가되는 시도 이력 |
| `collection_summary.json` | 선택 항목별 수집 결과와 catalog SHA256 |

- 중간 파일의 명세 fingerprint가 같을 때만 재개한다. URL/release/build/ancestry/checksum 등의 identity가 바뀌면 중단하고 새 대상 경로나 명세를 검토한다.
- `Range` 및 가능한 경우 `If-Range`로 재개한다. 206 응답의 Content-Range, 시작 offset, 길이와 ETag를 확인한다. 200 응답으로 Range를 무시하면 append하지 않고 처음부터 다시 쓴다.
- 짧게 끝난 HTTP body와 일시적 네트워크/408/429/5xx 오류는 지정된 횟수만 재시도한다. 중간 바이트는 남겨서 다음 실행도 재개할 수 있다.
- 최종 SHA256, 비어 있지 않은 파일, 선택적 크기를 확인한다. `.gz`/`.bgz`는 gzip 전체 stream 무결성과 비어 있지 않은 payload도 확인한다. checksum/압축 오류의 중간 파일은 폐기하고 실패를 기록한다.
- 기존 최종 파일은 덮어쓰지 않는다. 실제 바이트를 매번 재검증한다. 동일 config/code의 이전 검증 기록이면 `SKIPPED_EXISTING_VALID`, 이전 기록이 없거나 변경됐으면 `VERIFIED_EXISTING`으로 새 검증을 기록한다.
- 검증된 partial을 같은 디렉터리에 hard link로 최종 승격하므로 경합 시 이미 생긴 파일도 덮어쓰지 않는다. 파일과 메타데이터를 한 번에 묶는 transaction은 아니며, 승격 직후 프로세스가 끊겨도 다음 실행에서 최종 파일을 재검증한다.
- collection과 artifact에는 독점 lock이 있다. 정상 종료 시 제거한다. 프로세스가 강제 종료되어 lock이 남으면 실행 프로세스가 없다는 것을 확인한 후 해당 lock만 수동 처리한다. 자동으로 stale lock을 삭제하지 않는다.

`SUCCESS_VERIFIED`, `VERIFIED_EXISTING`, `SKIPPED_EXISTING_VALID`는 **파일 identity/integrity 상태**다. header-only 자료, allele/N/schema, genomic coordinates, study overlap, LD 적합성, scientific eligibility를 통과했다는 뜻이 아니다. 이 상태는 항상 `scientific_status=NOT_REVIEWED`와 함께 기록한다.

## 종료 코드와 기존 중앙 엔진 연결

- **0**: 선택된 모든 artifact의 integrity 검증 완료.
- **1**: 전송·로컬 파일·검증·잠금 등 실행 실패. batch에서는 항목별 실패를 기록하고 나머지 선택 항목은 처리한다.
- **2**: 하나 이상의 항목이 preflight에서 보류됨. batch 전체를 수집하지 않으며 대상 root도 만들지 않는다.

기존 `python -m masteromics acquire <dataset-spec.json> <marker.json>`도 원격 전송에 같은 `stage_file` 코어를 쓴다. 이미 있는 로컬 분석 입력은 같은 integrity 검증 함수를 호출하며 read-only 위치를 변경하지 않는다. 이 낮은 단계의 명세는 기존에 검토된 dataset registry를 사용하며, 신규 DB를 고를 때는 위 catalog/plan/collect 경로로 접근 조건부터 확인한다.

수집 summary의 `path`, `sha256`, release/build/ancestry를 후속 dataset spec의 provenance로 사용할 수 있다. column mapping과 형질/assay ID 선택은 아직 자동 생성하지 않는다. 현재 resource binding은 소스 ID당 한 파일이다. 여러 protein/assay/조직 파일은 검토된 artifact ID를 분리하는 registry 확장이 필요하다. S3 CLI, Synapse 로그인, tar member 추출 등의 서비스 전용 adapter는 이번 범위에 포함되지 않는다.

## CI

`database-catalog` job에서 NumPy/R 없이 기존 목록 검사와 수집기 테스트를 실행한다. 작은 `tests/fixtures/gigastroke.tsv`를 사용하며 외부 DB·credential은 필요 없다. 모의 응답과 실제 loopback HTTP server로 중단/재개를 확인한다. local controlled gate, checksum/gzip/range 실패, identity 변경, 기존 파일 보호, lock/symlink, 이력, stdlib CLI를 검사한다. 기존 Python/R 분석 CI도 함께 유지한다.
