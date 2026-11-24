#!/usr/bin/env python3
"""Source-level 18-library alignment, no expression counts or donor IDs exported."""
import argparse,collections,concurrent.futures,csv,gzip,hashlib,json,re,urllib.parse,urllib.request
from pathlib import Path
from ocular_gse313736_pool_pilot import META_SHA,norm,barcode_index,digest,rows,require

SERIES="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE313nnn/GSE313736"
HOST="ftp.ncbi.nlm.nih.gov"
ARMS=("no-dex","dex")
def fetch(url,path,cap):
    u=urllib.parse.urlsplit(url)
    require(u.scheme=="https" and u.hostname==HOST,"Unexpected download host")
    req=urllib.request.Request(url,headers={"User-Agent":"Public-ocular-expression-research/1.0"})
    with urllib.request.urlopen(req,timeout=40) as f,path.open("wb") as w:
        n=0
        while b:=f.read(1024*1024):
            n+=len(b)
            require(n<=cap,"Source download size bound")
            w.write(b)
    require(n>0,"Empty source")
    return digest(path)

def sample_records(path):
    records={};cur=None
    with gzip.open(path,"rt",encoding="utf-8-sig",errors="strict") as f:
        for line in f:
            line=line.rstrip("\r\n")
            if line.startswith("^"):
                cur=None
                if line.startswith("^SAMPLE = "):
                    cur=line.split(" = ",1)[1];require(re.fullmatch(r"GSM\d+",cur) and cur not in records,"Duplicate sample")
                    records[cur]=collections.defaultdict(list)
            elif cur and line.startswith("!Sample_") and " = " in line:
                key,value=line.split(" = ",1)
                records[cur][key].append(value)
    expected={f"GSM{n}" for n in range(9373276,9373294)}
    require(set(records)==expected,"Unexpected 18-library sample set")
    out={}
    for gsm,r in records.items():
        require("GSE313736" in r["!Sample_series_id"] and "GPL24676" in r["!Sample_platform_id"],"Sample/series/platform identity")
        desc=" ".join(r["!Sample_description"])
        m=re.search(r"Library name:\s*pool-(dex)?([A-I])\b",desc,re.I)
        require(m is not None,"Missing original pool name")
        arm="dex" if m.group(1) else "no-dex";pool=m.group(2).upper()
        expected_gsm=f"GSM{9373276+ord(pool)-ord('A')+(9 if arm=='dex' else 0)}"
        require(gsm==expected_gsm,"GSM/pool order differs")
        donors=set(map(int,re.findall(r"HTMC-(\d+)",desc,re.I)))
        require(bool(donors),"Missing original donor membership")
        excluded=set(map(int,re.findall(r"HTMC-(\d+)\s+was excluded",desc,re.I)))
        donors-=excluded
        stem="pool_"+("dex" if arm=="dex" else "")+pool
        files={}
        for url in r["!Sample_supplementary_file"]:
            url=url.replace("ftp://","https://",1)
            name=urllib.parse.urlsplit(url).path.rsplit("/",1)[-1]
            for kind,suffix in (("barcodes","barcodes.tsv.gz"),("features","features.tsv.gz"),("matrix","matrix.mtx.gz")):
                if name==gsm+"_"+stem+"_"+suffix:
                    require(kind not in files,"Repeated original supplemental path")
                    files[kind]=url
        require(set(files)=={"barcodes","features","matrix"},"Original supplementary file schema")
        duration_text=" ".join(r["!Sample_treatment_protocol"]+r["!Sample_growth_protocol"]+r["!Sample_characteristics_ch1"])
        durations=sorted(set(re.findall(r"\b(?:48\s*(?:hrs?|hours?)|3\s*days?)\b",duration_text,re.I)))
        out[pool,arm]={"gsm":gsm,"library":norm("pool-"+("dex" if arm=="dex" else "")+pool),
                       "donors":donors,"excluded":len(excluded),"files":files,"reported_duration_mentions":durations}
    require(len(out)==18,"Duplicate pool-arm records")
    return out

def main():
    p=argparse.ArgumentParser();p.add_argument("--root",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    root=a.root
    soft=root/"series.soft.gz";meta=root/"metadata.gz"
    soft_sha=fetch(SERIES+"/soft/GSE313736_family.soft.gz",soft,25000000)
    fetch(SERIES+"/suppl/GSE313736_HTMC_scRNA_cell_metadata.tsv.gz",meta,10000000)
    require(meta.stat().st_size==6592040 and digest(meta)==META_SHA,"Metadata pin")
    specs=sample_records(soft)
    labels={x["library"]:k for k,x in specs.items()}
    records=collections.defaultdict(list);locations=collections.defaultdict(dict);seen=set()
    with gzip.open(meta,"rt",encoding="utf-8-sig",newline="") as f:
        for r in csv.DictReader(f,delimiter="\t"):
            key=labels.get(norm(r["library_id"]))
            require(key is not None,"Unknown technical library label: "+r["library_id"])
            pool,arm=key;spec=specs[key]
            m=re.fullmatch(r"HTMC[-_]?(\d+)",r["sampleid"],re.I)
            require(m is not None and int(m.group(1)) in spec["donors"],"Source donor not in official pool")
            require(r["DEX Treatment"]==arm and r["droplet_type"]=="SNG","Arm/droplet mismatch")
            k=(key,r["barcode"]);require(k not in seen,"Duplicate library/barcode")
            seen.add(k);sid=r["sampleid"]
            require(arm not in locations[sid] or locations[sid][arm]==pool,"Sample in multiple pools per arm")
            locations[sid][arm]=pool
            records[key].append((sid,r["barcode"]))
    require(len(seen)==169398 and len(locations)==59 and set(records)==set(specs),"Fixed public cohort coverage")
    require(all(set(d)==set(ARMS) and d["dex"]==d["no-dex"] for d in locations.values()),"Cross-arm pool assignment differs")
    def verify(item):
        (pool,arm),spec=item
        b=root/f"{pool}-{arm}-barcodes.gz";f=root/f"{pool}-{arm}-features.gz"
        try:
            bs=fetch(spec["files"]["barcodes"],b,3000000);fs=fetch(spec["files"]["features"],f,3000000)
            bcs=rows(b);features=rows(f)
            require(all(len(x)==1 and re.fullmatch("[ACGT]+-[0-9]+",x[0]) for x in bcs),"Barcode schema")
            require(all(len(x)==3 and x[2]=="Gene Expression" and re.fullmatch(r"ENSG\d+(?:\.\d+)?",x[0]) for x in features),"Feature schema")
            ids=[x[0].split(".")[0] for x in features]
            require(len(set(ids))==len(ids),"Duplicate stable gene IDs")
            subset=records[pool,arm]
            mode,index=barcode_index([x[0] for x in bcs],[x[1] for x in subset])
            feature_hash=hashlib.sha256(("\n".join(ids)+"\n").encode()).hexdigest()
            return {"status":"EXACT_METADATA_BARCODE_FEATURE_ALIGNMENT","pool":pool,"arm":arm,"gsm":spec["gsm"],
                    "source_library_donors":len(spec["donors"]),"explicit_source_exclusions":spec["excluded"],
                    "metadata_samples":len({x[0] for x in subset}),"metadata_cells":len(subset),"matrix_barcodes":len(bcs),
                    "matched":len(index),"other_matrix_cells":len(bcs)-len(index),"barcode_mode":mode,
                    "feature_count":len(ids),"ordered_stable_feature_sha256":feature_hash,
                    "barcode_source_sha256":bs,"feature_source_sha256":fs,
                    "matrix_url":spec["files"]["matrix"],"reported_duration_mentions":spec["reported_duration_mentions"]}
        except (OSError,ValueError,EOFError,UnicodeError,csv.Error) as e:
            return {"status":"HOLD","pool":pool,"arm":arm,"gsm":spec["gsm"],"reason":str(e)[:300]}
        finally:b.unlink(missing_ok=True);f.unlink(missing_ok=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
        checks=list(ex.map(verify,sorted(specs.items())))
    pools=[]
    for pool in "ABCDEFGHI":
        d={arm:collections.Counter(s for s,b in records[pool,arm]) for arm in ARMS}
        require(set(d["dex"])==set(d["no-dex"]),"Unpaired metadata sample set")
        eligible={cut:sum(min(d[arm][s] for arm in ARMS)>=cut for s in d["dex"]) for cut in (50,100,200)}
        pools.append({"pool":pool,"metadata_pairs":len(d["dex"]),"pairs_by_min_cells":eligible,
                      "metadata_cells_by_arm":{arm:sum(d[arm].values()) for arm in ARMS},
                      "source_donor_set_equal":specs[pool,"dex"]["donors"]==specs[pool,"no-dex"]["donors"]})
    passed=[r for r in checks if r["status"]=="EXACT_METADATA_BARCODE_FEATURE_ALIGNMENT"]
    result={"status":"PUBLIC_18_LIBRARY_ALIGNMENT_GATE_NOT_EXPRESSION","soft_sha256":soft_sha,"metadata_sha256":META_SHA,
            "metadata_rows":len(seen),"unique_source_sample_ids":len(locations),"pools":pools,"libraries":checks,
            "passed":len(passed),"common_ordered_feature_universe":len(passed)==18 and len({r["ordered_stable_feature_sha256"] for r in passed})==1,
            "limits":["Barcode and feature alignment only, no counts scored. No donor identities or clinical data exported.",
                      "Matched released SNG metadata subset, not all published 74 donors/316799 cells.",
                      "Same-pool exposure; use pool-level inference rather than treating donors as independently treated replicates.",
                      "Source treatment/growth durations may disagree; verify original methods before reporting exact duration."]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps(result,indent=2,allow_nan=False))
    require(len(passed)==18 and result["common_ordered_feature_universe"],"Full library alignment gate failed")
if __name__=="__main__":main()
