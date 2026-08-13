#!/usr/bin/env python3
"""Phase 7B locked-TWPS spatial analysis. Plotting is deliberately delegated to R."""
from __future__ import annotations

import csv, hashlib, itertools, json, math, os, platform, shutil, sys, zipfile
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from scipy import sparse, stats

ROOT = Path(__file__).resolve().parents[1]
SP = ROOT / "06_results/spatial"
LOG = ROOT / "05_logs"
WORK = ROOT / "work/phase7b"
STAGED = ROOT / "data/staged/GSE241124/extracted"
MAP = SP / "GSE241124_SPATIAL_SAMPLE_MAP.tsv"
EXPECTED = {
 "CORE100": "6d5a62c080f2d0d9e4077983294d2a6a443a9e0e3c565691a4d3a863cc2afcc3",
 "CORE25": "fc9732278f494593640ba9088e0673a223648ba8f04d66b3c28d915051ee2ba0",
 "CORE50": "64717b1673393d2df74979843eaface382b39608975d7ac6414b43c74627ab2c",
 "FULL": "0110890f575ae60dbd177b610137a3c492057e4840d5586a0c42da83b7e649e1",
}
SIGFILES = {
 "CORE100": ROOT/"06_results/gateA/TWPS_PRIMARY_D7_M3_CORE100.tsv",
 "CORE25": ROOT/"06_results/gateA/TWPS_SENSITIVITY_D7_M3_CORE25.tsv",
 "CORE50": ROOT/"06_results/gateA/TWPS_SENSITIVITY_D7_M3_CORE50.tsv",
 "FULL": ROOT/"06_results/gateA/TWPS_SENSITIVITY_D7_M3_FULL.tsv",
}
TIME_ORDER = ["SKIN", "D1", "D7", "D30"]
CONTRASTS = [("D7_VS_SKIN","D7","SKIN",1),("D30_VS_D7","D30","D7",-1),
             ("D1_VS_SKIN","D1","SKIN",None),("D30_VS_SKIN","D30","SKIN",None)]
FIBRO = ["DCN","LUM","PDGFRA","COL1A1","COL1A2","COL3A1","COL6A1","COL6A2","COL6A3","DPT","CFD","PI16","FBLN1","FBLN2","COL15A1"]

def sha(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(1<<20),b""): h.update(b)
    return h.hexdigest()

def decode(a): return np.array([x.decode() if isinstance(x,bytes) else str(x) for x in a])

def read_h5(p):
    with h5py.File(p,"r") as h:
        g=h["matrix"]
        shape=tuple(int(x) for x in g["shape"][:])
        mat=sparse.csc_matrix((g["data"][:],g["indices"][:],g["indptr"][:]),shape=shape)
        bars=decode(g["barcodes"][:])
        names=decode(g["features/name"][:])
    return mat,names,bars

def read_spatial(zp, sample):
    dest=WORK/"spatial"/sample; dest.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(zp) as z:
        posname=next(n for n in z.namelist() if n.endswith("tissue_positions_list.csv"))
        sfname=next(n for n in z.namelist() if n.endswith("scalefactors_json.json"))
        imgname=next(n for n in z.namelist() if n.endswith("tissue_hires_image.png"))
        with z.open(posname) as f:
            pos=pd.read_csv(f,header=None,names=["barcode","in_tissue","array_row","array_col","pxl_row","pxl_col"])
        with z.open(sfname) as f: scales=json.load(f)
        img=dest/"tissue_hires_image.png"
        with z.open(imgname) as src, open(img,"wb") as out: shutil.copyfileobj(src,out)
    return pos,scales,img

def collapse_genes(mat,names):
    unq, inv=np.unique(names,return_inverse=True)
    if len(unq)==len(names): return mat,unq
    agg=sparse.csr_matrix((np.ones(len(inv)),(inv,np.arange(len(inv)))),shape=(len(unq),len(inv)))
    return agg@mat,unq

def exact_signflip(d):
    obs=abs(float(np.mean(d)))
    vals=[abs(float(np.mean(np.asarray(s)*d))) for s in itertools.product([-1,1],repeat=len(d))]
    return sum(v>=obs-1e-12 for v in vals)/len(vals)

def summarize(d, expected):
    d=np.asarray(d,float); n=len(d); mean=float(np.mean(d)); med=float(np.median(d)); sd=float(np.std(d,ddof=1))
    se=sd/math.sqrt(n); ci=stats.t.interval(.95,n-1,loc=mean,scale=se) if sd>0 else (mean,mean)
    dz=mean/sd if sd>0 else np.nan
    concord=int(np.sum(d>0)) if expected==1 else int(np.sum(d<0)) if expected==-1 else np.nan
    return dict(donor_n=n,mean_difference=mean,median_difference=med,sd=sd,ci_low=ci[0],ci_high=ci[1],effect_size=dz,pvalue=exact_signflip(d),direction_concordance=concord)

def main():
    SP.mkdir(parents=True,exist_ok=True); LOG.mkdir(exist_ok=True); WORK.mkdir(parents=True,exist_ok=True)
    hashes={k:sha(v) for k,v in SIGFILES.items()}
    if hashes != EXPECTED: raise SystemExit("TWPS_LOCK_MISMATCH")
    sigs={k:pd.read_csv(v,sep="\t")["gene"].astype(str).tolist() for k,v in SIGFILES.items()}
    core=set(sigs["CORE100"])
    smap=pd.read_csv(MAP,sep="\t")
    smap[["sample_id","GSM","donor_id","timepoint_standardized","included"]].to_csv(
        SP/"GSE241124_SLIDE_PSEUDOBULK_SAMPLE_MAP.tsv",sep="\t",index=False)
    rows=[]; spot_parts=[]; pseudo=[]; gene_union=None; slide_counts=[]; extracted=[]
    for r in smap.itertuples(index=False):
        h5name,zname=r.processed_file.split(";")
        hp,zp=STAGED/h5name,STAGED/zname
        mat,genes,bars=read_h5(hp); mat,genes=collapse_genes(mat,genes)
        pos,scales,img=read_spatial(zp,r.sample_id)
        pidx=pos.set_index("barcode"); keep=np.array([b in pidx.index and int(pidx.loc[b,"in_tissue"])==1 for b in bars])
        if not np.any(keep): raise RuntimeError(f"No in-tissue spots: {r.sample_id}")
        m=mat[:,keep].tocsc(); kb=bars[keep]; kp=pidx.loc[kb].reset_index()
        lib=np.asarray(m.sum(axis=0)).ravel(); detected=np.asarray((m>0).sum(axis=0)).ravel()
        counts=np.asarray(m.sum(axis=1)).ravel().astype(np.float64)
        slide_counts.append(pd.Series(counts,index=genes,name=r.sample_id))
        total_positions=len(pos); n_in=int(keep.sum())
        rows.append(dict(sample=r.sample_id,GSM=r.GSM,donor=r.donor_id,timepoint=r.timepoint_standardized,
                         h5_file=h5name,spatial_zip=zname,total_positions=total_positions,in_tissue_spots=n_in,
                         median_UMI=float(np.median(lib)),median_detected_genes=float(np.median(detected)),
                         library_size=float(lib.sum()),coordinate_barcode_match="YES",histology_image=str(img.relative_to(ROOT)),
                         h5_sha256=sha(hp),zip_sha256=sha(zp)))
        # library-size normalization and log1p on CPM for transparent spot scores
        norm=m.multiply(1e6/np.maximum(lib,1)).tocsc(); norm.data=np.log1p(norm.data)
        gidx={g:i for i,g in enumerate(genes)}
        def meanscore(glist):
            idx=[gidx[g] for g in glist if g in gidx]
            return np.asarray(norm[idx,:].mean(axis=0)).ravel(),len(idx)
        cscore,cavail=meanscore(sigs["CORE100"])
        z=(cscore-cscore.mean())/cscore.std(ddof=1) if cscore.std(ddof=1)>0 else np.zeros_like(cscore)
        independent=[g for g in FIBRO if g not in core]
        fscore,favail=meanscore(independent)
        q=np.quantile(fscore,.75); high=fscore>=q
        part=pd.DataFrame({"sample":r.sample_id,"donor":r.donor_id,"timepoint":r.timepoint_standardized,"barcode":kb,
                           "x":kp["pxl_col"].to_numpy(),"y":kp["pxl_row"].to_numpy(),"in_tissue":1,
                           "CORE100_spot_score":cscore,"SPOT_TWPS_DISPLAY_Z":z,"fibroblast_marker_score":fscore,
                           "fibroblast_enriched_top25":high.astype(int)})
        spot_parts.append(part)
        extracted.append(dict(sample=r.sample_id,hires_scale=float(scales["tissue_hires_scalef"]),histology_image=str(img.relative_to(ROOT))))
    qc=pd.DataFrame(rows); qc.to_csv(SP/"GSE241124_SPOT_QC.tsv",sep="\t",index=False)
    pd.DataFrame(rows).to_csv(LOG/"GSE241124_phase7B_manifest.tsv",sep="\t",index=False)
    pd.DataFrame(extracted).to_csv(WORK/"histology_manifest.tsv",sep="\t",index=False)
    # align all genes and make raw slide-pseudobulk matrix
    countdf=pd.concat(slide_counts,axis=1).fillna(0)
    countdf.index.name="gene"; countdf.to_csv(ROOT/"data/intermediate/GSE241124_slide_pseudobulk_counts.tsv.gz",sep="\t",compression="gzip")
    libs=countdf.sum(axis=0); logcpm=np.log2(countdf.divide(libs,axis=1)*1e6+1)
    logcpm.index.name="gene"; logcpm.to_csv(ROOT/"data/intermediate/GSE241124_slide_log2CPM.tsv.gz",sep="\t",compression="gzip")
    score_rows=[]; coverage=[]
    for sname,glist in sigs.items():
        avail=[g for g in glist if g in logcpm.index]
        zmat=logcpm.loc[avail].sub(logcpm.loc[avail].mean(axis=1),axis=0)
        sds=logcpm.loc[avail].std(axis=1,ddof=1); variable=sds[sds>0].index.tolist()
        zmat=zmat.loc[variable].div(sds.loc[variable],axis=0)
        scores=zmat.mean(axis=0)
        coverage.append(dict(signature=sname,locked_gene_n=len(glist),available_gene_n=len(avail),coverage=len(avail)/len(glist),variable_gene_n=len(variable)))
        for sample,val in scores.items():
            rr=smap.loc[smap.sample_id==sample].iloc[0]
            score_rows.append(dict(sample=sample,donor=rr.donor_id,timepoint=rr.timepoint_standardized,signature=sname,TWPS_score=val))
    scores=pd.DataFrame(score_rows); scores.to_csv(SP/"GSE241124_SLIDE_PSEUDOBULK_TWPS.tsv",sep="\t",index=False)
    pd.DataFrame(coverage).to_csv(SP/"GSE241124_SIGNATURE_COVERAGE_PHASE7B.tsv",sep="\t",index=False)
    scores[scores.signature=="CORE100"].sort_values(["donor","timepoint"]).to_csv(SP/"GSE241124_DONOR_TRAJECTORIES.tsv",sep="\t",index=False)
    summaries=[]; diffs=[]
    for sname in sigs:
        ss=scores[scores.signature==sname].pivot(index="donor",columns="timepoint",values="TWPS_score")
        for cname,a,b,expected in CONTRASTS:
            d=ss[a]-ss[b]; sm=summarize(d,expected)
            observed="POSITIVE" if sm["mean_difference"]>0 else "NEGATIVE" if sm["mean_difference"]<0 else "ZERO"
            exp="POSITIVE" if expected==1 else "NEGATIVE" if expected==-1 else "DESCRIPTIVE"
            summaries.append(dict(contrast=cname,signature=sname,**sm,expected_direction=exp,observed_direction=observed,evidence="PRESPECIFIED_PAIRED_DONOR_ANALYSIS" if expected else "SECONDARY_DESCRIPTIVE"))
            for donor,val in d.items(): diffs.append(dict(contrast=cname,signature=sname,donor=donor,difference=val))
    summ=pd.DataFrame(summaries); summ.to_csv(SP/"GSE241124_SPATIAL_TWPS_SUMMARY.tsv",sep="\t",index=False)
    pd.DataFrame(diffs).to_csv(SP/"GSE241124_PAIRED_DONOR_DIFFERENCES.tsv",sep="\t",index=False)
    summ[summ.signature!="CORE100"].to_csv(SP/"GSE241124_SIGNATURE_SENSITIVITY.tsv",sep="\t",index=False)
    nondisc=scores[(scores.signature=="CORE100") & scores.donor.isin(["Donor1","Donor2"])]
    ndwide=nondisc.pivot(index="donor",columns="timepoint",values="TWPS_score")
    nd=[]
    for donor,row in ndwide.iterrows():
        nd.append(dict(donor=donor,SKIN=row.SKIN,D1=row.D1,D7=row.D7,D30=row.D30,D7_GT_SKIN=row.D7>row.SKIN,D30_LT_D7=row.D30<row.D7))
    pd.DataFrame(nd).to_csv(SP/"GSE241124_NON_DISCOVERY_DONOR_SENSITIVITY.tsv",sep="\t",index=False)
    spots=pd.concat(spot_parts,ignore_index=True)
    spots[["sample","donor","timepoint","barcode","x","y","in_tissue","CORE100_spot_score","SPOT_TWPS_DISPLAY_Z"]].to_csv(SP/"GSE241124_SPOT_TWPS.tsv",sep="\t",index=False)
    retained=[g for g in FIBRO if g not in core]
    pd.DataFrame({"gene":FIBRO,"in_CORE100":[g in core for g in FIBRO],"retained_for_fibroblast_score":[g in retained for g in FIBRO]}).to_csv(SP/"GSE241124_FIBROBLAST_MARKER_PANEL.tsv",sep="\t",index=False)
    spots[["sample","barcode","fibroblast_marker_score"]].to_csv(SP/"GSE241124_FIBROBLAST_ENRICHMENT.tsv",sep="\t",index=False)
    assoc=[]
    for sample,g in spots.groupby("sample",sort=False):
        rho=stats.spearmanr(g.CORE100_spot_score,g.fibroblast_marker_score).statistic
        hi=g.fibroblast_enriched_top25.astype(bool)
        assoc.append(dict(sample=sample,donor=g.donor.iloc[0],timepoint=g.timepoint.iloc[0],spot_n=len(g),spearman_rho=rho,
                          top25_threshold=float(g.fibroblast_marker_score.quantile(.75)),top25_spot_n=int(hi.sum()),
                          top25_mean_TWPS=float(g.loc[hi,"CORE100_spot_score"].mean()),remainder_mean_TWPS=float(g.loc[~hi,"CORE100_spot_score"].mean()),
                          top25_minus_remainder_mean_TWPS=float(g.loc[hi,"CORE100_spot_score"].mean()-g.loc[~hi,"CORE100_spot_score"].mean())))
    pd.DataFrame(assoc).to_csv(SP/"GSE241124_TWPS_FIBROBLAST_ASSOCIATION.tsv",sep="\t",index=False)
    env={"python":sys.version.split()[0],"platform":platform.platform(),"numpy":np.__version__,"pandas":pd.__version__}
    (LOG/"GSE241124_phase7B_environment.json").write_text(json.dumps(env,indent=2)+"\n")
    print(summ.to_string(index=False)); print("FIBRO_RETAINED",len(retained)); print("SPOTS",len(spots))

if __name__=="__main__": main()
