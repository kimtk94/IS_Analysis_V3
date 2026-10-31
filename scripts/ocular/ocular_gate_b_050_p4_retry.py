#!/usr/bin/env python3
"""Gate B score .50 wrapper for P4 with bounded official GEO retry.

Does not change the scientific pseudobulk pipeline. Retries only transport
errors / invalid gzip before any count processing; original SHA256 and gzip
CRC, gene-order/barcode/NNZ checks remain enforced.
"""
import sys,time,urllib.parse,gzip
from ocular_p1_integrity import download as original
import ocular_p4_fullcohort as p4

def resilient_geo_download(url,path,cap):
    routes=[url]
    if "/suppl/" in url and "ftp.ncbi.nlm.nih.gov/geo/" in url:
        basename=url.rsplit("/",1)[-1]
        if basename.startswith("GSE"):
            accession=basename.split("_",1)[0]
        elif basename.startswith("GSM"):
            accession=basename.split("_",1)[0]
        else:
            accession=None
        if accession:
            routes.append("https://www.ncbi.nlm.nih.gov/geo/download/?"+
                urllib.parse.urlencode({"acc":accession,"format":"file","file":basename}))
    # Mirrors are attempted only after the canonical official source failed.
    last=None
    for route in routes:
        for attempt in range(3):
            try:
                metadata=original(route,path,cap)
                if metadata["bytes"]<=10:
                    raise ValueError("TRUNCATED_GEO_RESPONSE")
                if cap<=12*1024*1024:
                    # Validate complete gzip stream (not just magic) for small
                    # barcode/feature/cell metadata before downstream parsing.
                    with gzip.open(path,"rb") as src:
                        while src.read(1<<20):
                            pass
                if route!=url:
                    metadata["official_alternate_route_used"]=True
                return metadata
            except (OSError,ValueError,EOFError) as e:
                last=type(e).__name__+": "+str(e)[:170]
                print("GEO_RETRY",route.split("/")[2],"attempt",attempt+1,
                      "cause",last,flush=True)
                if attempt<2:time.sleep(2*(attempt+1))
    raise RuntimeError("ALL_OFFICIAL_GEO_DOWNLOAD_ROUTES_FAILED "+str(last))
def main():
    p4.download=resilient_geo_download
    p4.main()
if __name__=="__main__":main()
