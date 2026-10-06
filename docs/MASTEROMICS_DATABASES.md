# MasterOmics 데이터 소스와 접근 계획

검토일: 2026-10-03 UTC. 이 문서는 **자료 확보 경로와 실행 구조**를 정리한다. 실제 파일 다운로드, 헤더/행 검증, 수치 재현성, 독립 검증 코호트 확정은 이후 데이터 검토 단계에서 진행한다. 공개 포털이 있다는 사실과 현재 특정 파일을 받을 수 있다는 사실을 구분한다.

기준 목록은 `data/metadata/masteromics_resource_registry.json`이다. 분석 입력을 지정하는 `projects/datasets.example.json`과는 별개다. 목록의 resource/dataset ID는 소스 선택용이며, 포털 단위 항목은 정확한 연구·assay·조직·release·파일을 선택한 뒤 분석 입력 ID로 분리해야 한다.

## 우선 연결할 자료

| 자료군 | 확보 경로 | MasterOmics에서의 역할 | 접근/검토 조건 |
|---|---|---|---|
| CKDGen | [공식 데이터 목록](https://ckdgen.imbi.uni-freiburg.de/datasets/) | EUR eGFR/CKD, 보조 eGFRcys/BUN/UACR GWAS | 공개 요약통계. 기존 manifest의 5개 파일 경로를 보존하되 개별 endpoint는 이번에 다운로드하지 않았다. |
| TWB+BBJ kidney meta-GWAS | 기존 CKD manifest의 Figshare article 24356587 및 2개 파일 ID | EAS eGFR/BUN | 기존 경로를 가져온 후보. endpoint와 논문 메타데이터 재확인 필요. 포함된 BBJ를 독립 복제로 다시 세면 안 된다. |
| GIGASTROKE / MEGASTROKE | [GIGASTROKE 논문 Data availability](https://www.nature.com/articles/s41586-022-05165-3), [MEGASTROKE](https://www.megastroke.org/) | IS/하위유형 발견, ancestry별 분석, 과거 기준 비교 | GIGASTROKE GCST90104534–GCST90104563 중 정확한 형질/accession 선택. 두 연구의 코호트 중복 검토 필수. |
| UKB-PPP / deCODE | [UKB-PPP 공개 저장소 안내](https://registry.opendata.aws/ukbppp/), [deCODE summary data](https://www.decode.com/summarydata/) | pQTL 단백질 노출, MR·coloc·fine-mapping | 공개 요약통계. 플랫폼/assay, cis·전체 영역, ancestry를 선택. 개인 수준 단백질 자료 접근은 별도. |
| 1000 Genomes / IGSR | [공식 데이터](https://www.internationalgenome.org/data) | ancestry에 맞춘 reference LD | release/build/population/sample list 고정. EUR와 EAS를 분리하고 exposure/outcome LD 적합성 검토. |
| GTEx / eQTL Catalogue / NephQTL2 | [GTEx](https://gtexportal.org/home/downloads/adult-gtex/), [eQTL Catalogue](https://www.ebi.ac.uk/eqtl/Data_access/), [NephQTL2 배포 안내](https://www.sampsonlab.org/tools-and-data-sets) | 조직·세포별 조절 유전 근거 | 처리된 QTL은 공개 경로가 있다. 유의 SNP만 있는 파일은 full-region coloc 입력으로 부족하다. |
| HPA / Allen human brain / KPMP | [HPA](https://www.proteinatlas.org/about/download), [Allen WHB](https://alleninstitute.github.io/abc_atlas_access/descriptions/WHB-10Xv3.html), [KPMP](https://www.kpmp.org/available-data) | 단백질/조직/세포 위치 주석 | Allen 정상 사람 뇌, KPMP reference/AKI/CKD를 구분. 발현과 질환 인과성은 별도 근거다. |
| BBJ / FinnGen / GWAS Catalog | [BBJ](https://pheweb.jp/downloads), [FinnGen](https://www.finngen.fi/en/access_results), [GWAS Catalog](https://www.ebi.ac.uk/gwas/docs/file-downloads/) | 추가 outcome·sensitivity·검증 후보 | BBJ/Catalog는 연구별 파일 선택. FinnGen은 양식 제출 후 이메일 다운로드 안내. 독립성은 소스 이름만으로 판단하지 않는다. |
| KoGES / UKB / GTEx 개인 자료 | [KoGES](https://www.nih.go.kr/eng/main/contents.do?menuNo=200079), [UKB](https://www.ukbiobank.ac.uk/use-our-data/), [GTEx 보호 자료 안내가 있는 페이지](https://gtexportal.org/home/downloads/adult-gtex/) | score, longitudinal eGFR, incident outcome 등 | 승인·사용계약·허용 범위 확인 후 별도 adapter. 공개 요약통계만으로 개인 수준 분석을 실행할 수 없다. |

권장 연결 순서는 **pQTL + CKD/IS GWAS + LD → MR/coloc/SuSiE → 추가 QTL/조직 주석 → 독립성 확인을 거친 복제 → 승인된 개인 수준 분석**이다. DB 수를 늘리는 것보다 정확한 파일과 역할을 먼저 고정한다.

## 접근 분류와 기록 범위

- `open`: 공개 자료 경로가 확인된 소스. 모든 파일의 재배포/상업 이용 허용을 뜻하지 않는다.
- `registration`: 가입·폼·토큰 등 사용자 절차가 필요하다. FinnGen, OpenGWAS API, FinnGen multiome 후보를 포함한다.
- `controlled`: 승인된 연구/DUA 등의 절차가 필요한 개인 수준 자료. 목록의 binding으로 이 접근 분류를 바꿀 수 없다. 이미 확보한 로컬 파일은 권한 확인 reference를 명시한 경우에 한해 별도 staging할 수 있다.
- `publication_only`, `unknown`: 정확한 자료 또는 접근 경로가 검증되지 않은 후보에 사용할 상태다.

`last_verified`는 위 접근 안내 또는 기존 manifest를 검토한 날짜다. `verification_scope=official_portal_review`는 공식 안내 검토이며, `existing_manifest_not_probed`는 기존 경로를 보존했으나 개별 파일을 확인하지 않았다는 뜻이다. 모든 항목은 `download_status=not_checked` 또는 `restricted`, `checksum_status=not_checked`로 시작한다. 미확인 N/build/ancestry는 임의로 채우지 않는다. 발표된 전체 N과 실제 SNP별 N은 구별한다.

GTEx 공개 처리 자료와 보호된 원시/개인 자료는 별도 ID다. eQTL Catalogue는 공식 안내상 **REST API가 폐기되었으므로 FTP/파일과 release metadata**를 기준으로 adapter를 만든다. FinnGen 최신 release를 자동 선택하지 않는다. OpenGWAS JWT는 저장소·CI에 넣지 않는다.

사람 뇌졸중 환자 sc/snRNA, 동물 MCAO, spatial/ATAC 질환 자료는 현재 확정 accession이 없는 관계로 다운로드 가능하다고 등록하지 않았다. 추후 논문/GEO/SRA/공식 저장소에서 accession·donor·상태·조직·시점·종·공개 범위를 확인한 뒤 추가한다. 동물 결과를 사람 stroke differential expression으로 해석하지 않는다. KPMP의 개별 원시 자료 접근은 collection별로 다시 검토한다.

## 오프라인 CLI

NumPy/R 없이 목록과 계획을 검사할 수 있다. 서버에서 저장소 루트로 이동한 뒤 실행한다.

```bash
python3 -S -m masteromics resources validate
python3 -S -m masteromics resources list --project ckd
python3 -S -m masteromics resources list --project ischemic_stroke
python3 -S -m masteromics resources plan --ids ukb_ppp_pqtl gigastroke onekg_ld
```

`plan`은 다운로드하지 않는다. 미고정 항목이 있으면 JSON에 이유를 출력하고 exit **2**로 종료한다. 오류는 **1**, 모든 선택 항목의 접근/고정 조건을 충족하면 **0**이다. `READY_FOR_ACQUISITION`은 오프라인 메타데이터 조건 충족이며 endpoint/파일/생물학 검증 통과가 아니다.

정확한 공개 파일을 검토한 뒤에만 별도의 로컬 bindings JSON을 전달한다. 예시 필드는 다음과 같다. 아래 null을 실제 값으로 고정해야 하며 SHA256은 확인된 파일 또는 신뢰할 수 있는 출처의 digest를 사용한다.

```json
{
  "ukb_ppp_pqtl": {
    "artifact_url": null,
    "file_name": null,
    "release": null,
    "genome_build": null,
    "ancestry": null,
    "sha256": null,
    "license_reviewed": false
  }
}
```

```bash
python3 -S -m masteromics resources plan \
  --ids ukb_ppp_pqtl --bindings /absolute/path/resources.bindings.json
```

이 `plan` 명령은 acquisition을 실행하지 않는다. 검토된 binding으로 `resources collect`를 실행하면 공통 수집기로 파일과 이력을 저장할 수 있다. [수집기 사용법](MASTEROMICS_ACQUISITION.md)을 참고한다. 파일별 column adapter와 분석 dataset registry로 연결하는 단계는 남아 있다. 기존 CKD downloader는 `scripts/download_ckd_public_data.py`, UKB-PPP 파일 목록은 `data/metadata/ukb_ppp_download_manifest.tsv`에 있다. 공개 HTTPS 파일과 승인된 로컬 파일은 공통 수집기로 처리한다. 서비스 로그인/API 및 DB별 변환 adapter는 미구현이다.

## CI 검증

`database-catalog` job은 scientific 패키지 설치 없이 `python -S`로 실행한다. 실제 대용량 파일, 로그인/토큰, 외부 API 요청 없이 다음을 검사한다.

1. 필수 필드, 유일 ID, HTTPS URL, 안전한 파일명, 프로젝트·category·access 값.
2. 기존 CKD 공개 manifest의 파일 경로 보존과 검토 범위 표시.
3. exact artifact/release/SHA256/license/build/ancestry 미고정 시 acquisition plan 차단.
4. registration/controlled 항목의 binding으로 접근 분류를 우회할 수 없는지.
5. 사람 reference와 동물 category 혼동, 알 수 없는 ID, 중복 선택, 잘못된 binding 거부.
6. NumPy/Pandas 없는 환경에서도 CLI 동작 및 blocked plan exit 2.

기존 `python-contracts`, `r-syntax`, 실제 `r-integration` job도 유지한다. 실제 파일의 헤더·allele·N·checksum·LD·cohort overlap 검토와 서버 수치 비교는 별도 단계다. 외부 포털 장애를 오프라인 CI 성공/실패 판정에 섞지 않는다.

## 등록 ID 전체 목록

아래 표는 기준 JSON에서 생성했다. P1은 핵심 연결, P2는 추가 근거/검증 후보, P3는 이후 확장이다.

| ID | 소스 | 접근 | 프로젝트 | 우선순위 |
|---|---|---|---|---|
| `ckdgen_stanzick2021_egfr_eur` | [CKDGen eGFRcrea EUR](https://ckdgen.imbi.uni-freiburg.de/datasets/) | open | ckd | P1 |
| `ckdgen_wuttke2019_ckd_eur` | [CKDGen CKD EUR](https://ckdgen.imbi.uni-freiburg.de/datasets/) | open | ckd | P1 |
| `ckdgen_wuttke2019_bun_eur` | [CKDGen BUN EUR](https://ckdgen.imbi.uni-freiburg.de/datasets/) | open | ckd | P2 |
| `ckdgen_gorski2017_egfrcys_eur` | [CKDGen eGFRcys EUR](https://ckdgen.imbi.uni-freiburg.de/datasets/) | open | ckd | P2 |
| `ckdgen_teumer2019_uacr_eur` | [CKDGen UACR EUR](https://ckdgen.imbi.uni-freiburg.de/datasets/) | open | ckd | P3 |
| `eas_chen2024_egfr_meta` | [Figshare eGFRcrea EAS](https://figshare.com/articles/dataset/24356587) | open | ckd | P1 |
| `eas_chen2024_bun_meta` | [Figshare BUN EAS](https://figshare.com/articles/dataset/24356587) | open | ckd | P2 |
| `gigastroke` | [GIGASTROKE GWAS](https://www.nature.com/articles/s41586-022-05165-3) | open | ischemic_stroke | P1 |
| `megastroke` | [MEGASTROKE GWAS](https://www.megastroke.org/) | open | ischemic_stroke | P1 |
| `bbj_gwas` | [Biobank Japan / PheWeb](https://pheweb.jp/downloads) | open | ckd, ischemic_stroke | P1 |
| `finngen_gwas` | [FinnGen summary statistics](https://www.finngen.fi/en/access_results) | registration | ckd, ischemic_stroke | P2 |
| `gwas_catalog` | [NHGRI-EBI GWAS Catalog summary statistics](https://www.ebi.ac.uk/gwas/docs/file-downloads/) | open | ckd, ischemic_stroke | P2 |
| `ukb_ppp_pqtl` | [UKB-PPP public pGWAS](https://registry.opendata.aws/ukbppp/) | open | ckd, ischemic_stroke | P1 |
| `decode_pqtl` | [deCODE summary data](https://www.decode.com/summarydata/) | open | ckd, ischemic_stroke | P2 |
| `gtex_eqtl` | [GTEx adult processed QTL](https://gtexportal.org/home/downloads/adult-gtex/) | open | ckd, ischemic_stroke | P2 |
| `eqtl_catalogue` | [eQTL Catalogue](https://www.ebi.ac.uk/eqtl/Data_access/) | open | ckd, ischemic_stroke | P2 |
| `nephqtl2` | [NephQTL2](https://www.sampsonlab.org/tools-and-data-sets) | open | ckd | P2 |
| `onekg_ld` | [1000 Genomes / IGSR](https://www.internationalgenome.org/data) | open | ckd, ischemic_stroke | P1 |
| `hpa_annotation` | [Human Protein Atlas](https://www.proteinatlas.org/about/download) | open | ckd, ischemic_stroke | P2 |
| `allen_human_brain` | [Allen ABC Atlas human brain](https://alleninstitute.github.io/abc_atlas_access/descriptions/WHB-10Xv3.html) | open | ischemic_stroke | P2 |
| `kpmp_processed` | [KPMP processed kidney atlas](https://www.kpmp.org/available-data) | open | ckd | P2 |
| `koges_cohort` | [KoGES individual-level cohort](https://www.nih.go.kr/eng/main/contents.do?menuNo=200079) | controlled | ckd, ischemic_stroke | P3 |
| `ukb_individual` | [UK Biobank individual-level data](https://www.ukbiobank.ac.uk/use-our-data/) | controlled | ckd, ischemic_stroke | P3 |
| `gtex_controlled` | [GTEx individual-level sequencing/genotypes](https://gtexportal.org/home/downloads/adult-gtex/) | controlled | ckd, ischemic_stroke | P3 |
| `opengwas` | [IEU OpenGWAS API](https://api.opengwas.io/api/docs) | registration | ckd, ischemic_stroke | P3 |
| `finngen_multiome` | [FinnGen immune multiome QTL](https://www.finngen.fi/en/access_results) | registration | ckd, ischemic_stroke | P3 |
