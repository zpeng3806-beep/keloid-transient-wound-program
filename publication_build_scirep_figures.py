#!/usr/bin/env python3
"""Phase 10B presentation-only figure builder.

Reads frozen project outputs, performs no model fitting or biological inference,
and renders with one Python-only vector stack: ReportLab for SVG/PDF and Pillow
for high-resolution PNG diagnostics.
"""
from pathlib import Path
import argparse, csv, math, hashlib, re, shutil, copy
from reportlab.graphics.shapes import Drawing, Rect, Line, Circle, String, Polygon, PolyLine
from reportlab.graphics import renderSVG, renderPDF
from reportlab.lib.colors import HexColor, Color
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from PIL import Image, ImageDraw, ImageFont
from pypdf import PdfReader, PdfWriter
from pypdf.generic import DecodedStreamObject

ROOT = Path(__file__).resolve().parents[3]
PUB = ROOT / "08_publication"
OUT = Path(__file__).resolve().parent
MAIN = OUT / "figure_exports"
SUPP = OUT / "supplementary_figure_exports"
TABLES = PUB / "tables"
SOURCE = PUB / "source_data"
FONT_REGULAR_PATH = "/System/Library/Fonts/Supplemental/Arial.ttf"
FONT_BOLD_PATH = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"
FONT_REGULAR = "ArialEmbedded"
FONT_BOLD = "ArialEmbedded-Bold"
pdfmetrics.registerFont(TTFont(FONT_REGULAR, FONT_REGULAR_PATH))
pdfmetrics.registerFont(TTFont(FONT_BOLD, FONT_BOLD_PATH))

def remove_unused_initial_pdf_font(path):
    """Remove ReportLab's empty Times-Roman setup commands and resource."""
    path=Path(path)
    writer=PdfWriter(clone_from=path)
    for page in writer.pages:
        content=page.get_contents().get_data()
        content=re.sub(rb"BT\s*/F1\s+10\s+Tf\s+12\s+TL\s+ET\s*",b"",content)
        stream=DecodedStreamObject()
        stream.set_data(content)
        page.replace_contents(stream)
        fonts=page["/Resources"].get("/Font")
        if fonts:
            fonts=fonts.get_object()
        if fonts and b"/F1" not in content and "/F1" in fonts:
            del fonts["/F1"]
    temp=path.with_suffix(".fontfix.pdf")
    with open(temp,"wb") as handle:
        writer.write(handle)
    temp.replace(path)

# Vector renderer equivalents of matplotlib's required publication settings:
# svg.fonttype = 'none'; pdf.fonttype = 42; font.size = 7.
# FINAL_WIDTH_MM = 180; PNG_DPI = 600. Text remains native/editable in SVG/PDF.

PAL = {
    "Skin": "#5B6F82", "UNINJURED_SKIN": "#5B6F82", "D1": "#E69F00",
    "EARLY_WOUND": "#E69F00", "D7": "#8E63A9", "D30": "#2A9D8F",
    "LATE_WOUND": "#2A9D8F", "healthy": "#4C78A8", "keloid-prone": "#CC79A7",
    "normal_scar": "#7F8C8D", "keloid": "#6A3D9A", "HYPERTROPHIC_SCAR": "#D9A441",
    "Other": "#B8BDC5", "dark": "#28323C", "light": "#F4F6F8", "white": "#FFFFFF",
    "negative": "#3D6FA3", "positive": "#B55A5A",
}

def rows(path):
    with open(ROOT / path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def fnum(v):
    try: return float(v)
    except Exception: return math.nan

def subset(data, **conds):
    return [r for r in data if all(r.get(k) == v for k, v in conds.items())]

class Scene:
    def __init__(self, w=510, h=650, title=""):
        self.w, self.h = w, h
        self.d = Drawing(w, h)
        self.ops = []
        self.rect(0, 0, w, h, fill="#FFFFFF", stroke=None)
        if title:
            self.text(18, h-18, title, 11, bold=True, anchor="start")
            self.line(18, h-24, w-18, h-24, "#C8CDD2", 0.6)

    def rect(self, x,y,w,h,fill=None,stroke="#444444",sw=0.7,rx=0):
        self.d.add(Rect(x,y,w,h,rx=rx,ry=rx,fillColor=HexColor(fill) if fill else None,
                        strokeColor=HexColor(stroke) if stroke else None,strokeWidth=sw))
        self.ops.append(("rect",x,y,w,h,fill,stroke,sw))

    def line(self,x1,y1,x2,y2,color="#444444",sw=0.8,dash=None):
        obj=Line(x1,y1,x2,y2,strokeColor=HexColor(color),strokeWidth=sw)
        if dash: obj.strokeDashArray=dash
        self.d.add(obj); self.ops.append(("line",x1,y1,x2,y2,color,sw,dash))

    def circle(self,x,y,r,fill="#777777",stroke="#FFFFFF",sw=0.35):
        self.d.add(Circle(x,y,r,fillColor=HexColor(fill),strokeColor=HexColor(stroke) if stroke else None,strokeWidth=sw))
        self.ops.append(("circle",x,y,r,fill,stroke,sw))

    def polyline(self,pts,color="#444444",sw=1.0):
        flat=[z for p in pts for z in p]
        self.d.add(PolyLine(flat,strokeColor=HexColor(color),strokeWidth=sw,fillColor=None))
        self.ops.append(("polyline",pts,color,sw))

    def polygon(self,pts,fill="#DDDDDD",stroke="#444444",sw=0.7):
        flat=[z for p in pts for z in p]
        self.d.add(Polygon(flat,fillColor=HexColor(fill),strokeColor=HexColor(stroke),strokeWidth=sw))
        self.ops.append(("polygon",pts,fill,stroke,sw))

    def text(self,x,y,s,size=7,color="#28323C",bold=False,anchor="start",angle=0):
        size=max(7,size)
        font=FONT_BOLD if bold else FONT_REGULAR
        self.d.add(String(x,y,str(s),fontName=font,fontSize=size,fillColor=HexColor(color),textAnchor=anchor,angle=angle))
        self.ops.append(("text",x,y,str(s),size,color,bold,anchor,angle))

    def wrapped(self,x,y,s,width_chars=42,size=7,leading=9,color="#28323C",bold=False):
        words=str(s).split(); lines=[]; cur=""
        for word in words:
            test=(cur+" "+word).strip()
            if len(test)>width_chars and cur: lines.append(cur); cur=word
            else: cur=test
        if cur: lines.append(cur)
        for i,line in enumerate(lines): self.text(x,y-i*leading,line,size,color,bold)
        return y-len(lines)*leading

    def panel(self,x,y,w,h,label,title=""):
        self.rect(x,y,w,h,fill="#FFFFFF",stroke="#D5D9DD",sw=0.55,rx=2)
        self.text(x+6,y+h-12,label.lower(),8,bold=True)
        if title: self.text(x+20,y+h-12,title,7.2,bold=True)

    def export(self, stem):
        stem=Path(stem); stem.parent.mkdir(parents=True,exist_ok=True)
        renderSVG.drawToFile(self.d,str(stem.with_suffix('.svg')))
        pdf_path=stem.with_suffix('.pdf')
        renderPDF.drawToFile(self.d,str(pdf_path),initialFontName=FONT_REGULAR)
        remove_unused_initial_pdf_font(pdf_path)
        # ReportLab's PostScript driver cannot embed TrueType fonts. Use an
        # otherwise identical vector scene with native Helvetica fonts for EPS;
        # SVG/PDF retain embedded Arial. No raster image is introduced.
        from reportlab.graphics import renderPS
        eps_drawing = copy.deepcopy(self.d)
        def eps_fonts(node):
            for child in getattr(node, 'contents', []):
                if isinstance(child, String):
                    child.fontName = 'Helvetica-Bold' if 'Bold' in child.fontName else 'Helvetica'
                eps_fonts(child)
        eps_fonts(eps_drawing)
        renderPS.drawToFile(eps_drawing, str(stem.with_suffix('.eps')), fmt='EPS')
        scale=8
        im=Image.new("RGB",(int(self.w*scale),int(self.h*scale)),"white")
        dr=ImageDraw.Draw(im)
        def xy(x,y): return (int(x*scale),int((self.h-y)*scale))
        try: regular=ImageFont.truetype(FONT_REGULAR_PATH,7*scale)
        except Exception: regular=ImageFont.load_default()
        fonts={}
        for op in self.ops:
            if op[0]=="rect":
                _,x,y,w,h,fill,stroke,sw=op; a=xy(x,y+h); b=xy(x+w,y)
                dr.rectangle([a,b],fill=fill,outline=stroke,width=max(1,int(sw*scale)) if stroke else 0)
            elif op[0]=="line":
                _,x1,y1,x2,y2,c,sw,dash=op; dr.line([xy(x1,y1),xy(x2,y2)],fill=c,width=max(1,int(sw*scale)))
            elif op[0]=="circle":
                _,x,y,r,fill,stroke,sw=op; a=xy(x-r,y+r); b=xy(x+r,y-r)
                dr.ellipse([a,b],fill=fill,outline=stroke,width=max(1,int(sw*scale)) if stroke else 0)
            elif op[0]=="polyline":
                _,pts,c,sw=op; dr.line([xy(*p) for p in pts],fill=c,width=max(1,int(sw*scale)),joint="curve")
            elif op[0]=="polygon":
                _,pts,fill,stroke,sw=op; dr.polygon([xy(*p) for p in pts],fill=fill,outline=stroke)
            elif op[0]=="text":
                _,x,y,s,size,c,bold,anchor,angle=op
                key=(int(size*scale),bold)
                if key not in fonts:
                    path=FONT_BOLD_PATH if bold else FONT_REGULAR_PATH
                    try: fonts[key]=ImageFont.truetype(path,max(6,int(size*scale)))
                    except Exception: fonts[key]=regular
                font=fonts[key]; px,py=xy(x,y); box=dr.textbbox((0,0),s,font=font); tw=box[2]
                if anchor=="middle": px-=tw//2
                elif anchor=="end": px-=tw
                dr.text((px,py-int(size*scale)),s,font=font,fill=c)
        im.save(stem.with_suffix('.png'),dpi=(600,600),optimize=True)

def axes(sc,x,y,w,h,ymin,ymax,xlabels=None,ylabel=""):
    sc.line(x,y,x+w,y,"#555555",0.7); sc.line(x,y,x,y+h,"#555555",0.7)
    if ymin<0<ymax:
        yy=y+(0-ymin)/(ymax-ymin)*h; sc.line(x,yy,x+w,yy,"#AEB4BA",0.5,[2,2])
    for k in range(3):
        v=ymin+(ymax-ymin)*k/2; yy=y+h*k/2
        sc.line(x-2,yy,x,yy,"#555555",0.5); sc.text(x-4,yy-2,f"{v:.1f}",5.6,anchor="end")
    if xlabels:
        pretty={'UNINJURED_SKIN':'Uninjured','EARLY_WOUND':'Early','LATE_WOUND':'Late','HYPERTROPHIC_SCAR':'HTS','normal_scar':'Normal scar','keloid-prone':'Keloid-prone'}
        for i,l in enumerate(xlabels):
            xx=x+w*(i+0.5)/len(xlabels); sc.text(xx,y-10,pretty.get(l,l.replace("_"," ").title()),5.8,anchor="middle")
    if ylabel: sc.text(x,y+h+5,ylabel,5.8,anchor="start")
    return lambda v: y+(v-ymin)/(ymax-ymin)*h

def dot_groups(sc,box,data,group_key,value_key,groups,colors=None,ymin=None,ymax=None,ylabel="TWPS"):
    x,y,w,h=box; vals=[fnum(r[value_key]) for r in data if r.get(group_key) in groups]
    ymin=min(vals) if ymin is None else ymin; ymax=max(vals) if ymax is None else ymax
    pad=max((ymax-ymin)*.12,.1); ymin-=pad; ymax+=pad
    yy=axes(sc,x,y,w,h,ymin,ymax,groups,ylabel)
    for gi,g in enumerate(groups):
        rs=[r for r in data if r.get(group_key)==g]
        cx=x+w*(gi+.5)/len(groups)
        for j,r in enumerate(rs):
            jitter=((int(hashlib.sha256((g+str(j)).encode()).hexdigest()[:4],16)%101)/100-0.5)*min(18,w/len(groups)*.35)
            sc.circle(cx+jitter,yy(fnum(r[value_key])),2.0,PAL.get(g,PAL["Other"]),"#FFFFFF",.25)

def paired_lines(sc,box,data,id_key,time_key,value_key,times,colors=None,ylabel="TWPS"):
    x,y,w,h=box; vals=[fnum(r[value_key]) for r in data if r.get(time_key) in times]
    pad=max((max(vals)-min(vals))*.12,.1); ymin=min(vals)-pad; ymax=max(vals)+pad
    yy=axes(sc,x,y,w,h,ymin,ymax,times,ylabel)
    ids=[]
    for r in data:
        if r[id_key] not in ids: ids.append(r[id_key])
    for idx,ident in enumerate(ids):
        rs={r[time_key]:r for r in data if r[id_key]==ident}
        pts=[]
        for ti,t in enumerate(times):
            if t in rs:
                xx=x+w*(ti+.5)/len(times); v=yy(fnum(rs[t][value_key])); pts.append((xx,v))
        if len(pts)>1: sc.polyline(pts,"#87919A",0.65)
        for ti,(xx,v) in enumerate(pts): sc.circle(xx,v,2.1,(colors or PAL).get(times[ti],PAL["Other"]),"#FFFFFF",.25)

def bar_values(sc,box,labels,values,colors=None,ylabel="",zero=True):
    x,y,w,h=box; lo=min(values+[0]) if zero else min(values); hi=max(values+[0]); span=max(hi-lo,1e-6)
    if zero and lo>=0: lo=0
    else: lo-=span*.12
    if zero and hi<=0: hi=0
    else: hi+=span*.16
    yy=axes(sc,x,y,w,h,lo,hi,labels,ylabel)
    bw=w/len(labels)*.52; z=yy(0 if zero else lo)
    for i,(lab,v) in enumerate(zip(labels,values)):
        cx=x+w*(i+.5)/len(labels); yv=yy(v); col=(colors[i] if colors else (PAL["positive"] if v>=0 else PAL["negative"]))
        sc.rect(cx-bw/2,min(z,yv),bw,abs(yv-z),fill=col,stroke=None)
        sc.text(cx,yv+(3 if v>=0 else -7),f"{v:.2f}",5.3,anchor="middle")

def forest(sc,box,labels,est,low,high,colors=None,xlabel="Difference"):
    x,y,w,h=box; lo=min(low+[0]); hi=max(high+[0]); pad=(hi-lo)*.12; lo-=pad; hi+=pad
    xx=lambda v:x+(v-lo)/(hi-lo)*w
    sc.line(xx(0),y,xx(0),y+h,"#AEB4BA",0.6,[2,2])
    for i,lab in enumerate(labels):
        yy=y+h*(len(labels)-i-.5)/len(labels); sc.text(x-4,yy-2,lab,5.7,anchor="end")
        sc.line(xx(low[i]),yy,xx(high[i]),yy,"#4A5560",1.0); sc.line(xx(low[i]),yy-2,xx(low[i]),yy+2,"#4A5560",.6); sc.line(xx(high[i]),yy-2,xx(high[i]),yy+2,"#4A5560",.6)
        sc.circle(xx(est[i]),yy,2.7,(colors[i] if colors else PAL["dark"]),"#FFFFFF",.3)
    sc.line(x,y,x+w,y,"#555555",.6); sc.text(x+w/2,y-11,xlabel,5.8,anchor="middle")

def spatial(sc,box,data,sample,title):
    x,y,w,h=box; rs=[r for r in data if r["sample"]==sample and r["in_tissue"]=="1"]
    if not rs:
        sc.text(x+w/2,y+h/2,"No map data",7,anchor="middle"); return
    xs=[fnum(r["x"]) for r in rs]; ys=[fnum(r["y"]) for r in rs]
    xmin,xmax=min(xs),max(xs); ymin,ymax=min(ys),max(ys)
    def col(v):
        v=max(-2,min(2,v)); t=(v+2)/4
        a=(55,88,135); b=(230,159,0); rgb=tuple(int(a[i]*(1-t)+b[i]*t) for i in range(3)); return '#%02X%02X%02X'%rgb
    for r in rs:
        xx=x+(fnum(r['x'])-xmin)/(xmax-xmin)*w; yy=y+(fnum(r['y'])-ymin)/(ymax-ymin)*h
        sc.circle(xx,yy,.75,col(fnum(r['SPOT_TWPS_DISPLAY_Z'])),None,0)
    sc.text(x+w/2,y+h+4,title,5.8,bold=True,anchor="middle")

def figure1():
    sc=Scene()
    boxes=[(18,390,474,220),(18,205,150,165),(180,205,150,165),(342,205,150,165),(18,28,231,157),(261,28,231,157)]
    for b,l,t in zip(boxes,"ABCDEF",["Study design and lock","GSE241132 design","Transient discovery","Five modules","D7_M3 trajectory","CORE100 lock"]): sc.panel(*b,l,t)
    x,y,w,h=boxes[0]; stages=[("Normal wound\ndiscovery",PAL['D1']),("Candidate\nprograms",PAL['D7']),("D7_M3 +\nCORE100",PAL['D7']),("Disease and spatial\nvalidation",PAL['keloid']),("Prospective\ngenetics",PAL['Other'])]
    sx=x+14; cy=y+92
    for i,(lab,c) in enumerate(stages):
        bw=72 if i<3 else 84; sc.rect(sx,cy,bw,42,fill=c,stroke=None,rx=5); sc.wrapped(sx+bw/2,cy+26,lab,15,6,7,"#FFFFFF",True); sx+=bw+15
        if i<len(stages)-1: sc.line(sx-12,cy+21,sx-3,cy+21,"#616970",1.2)
    lockx=x+267; sc.rect(lockx,y+34,5,145,fill="#1F2D3D",stroke=None); sc.text(lockx+7,y+176,"TWPS LOCKED",7,bold=True); sc.text(lockx+7,y+166,"before disease outcomes",6)
    sc.text(x+14,y+42,"TWPS = Transient Wound Program Score",6.2,bold=True)
    sm=rows("06_results/gateA/PSEUDOBULK_SAMPLE_MAP.tsv"); x,y,w,h=boxes[1]
    for di,d in enumerate(sorted(set(r['donor'] for r in sm))):
        sc.text(x+12,y+h-35-di*34,d,6,bold=True)
        for ti,t in enumerate(['Skin','D1','D7','D30']): sc.circle(x+52+ti*23,y+h-33-di*34,5,PAL[t],None,0)
    sc.text(x+w/2,y+14,"3 donors × 4 timepoints",6.2,anchor="middle")
    summ={r['metric']:r['value'] for r in rows("06_results/gateA/gate_a2_summary.tsv")}; x,y,w,h=boxes[2]
    sc.text(x+15,y+h-42,"Tier 2",6); sc.text(x+w-15,y+h-42,summ['TIER2_N'],12,bold=True,anchor="end")
    sc.text(x+15,y+h-80,"Robust core",6); sc.text(x+w-15,y+h-80,summ['ROBUST_CORE_N'],12,bold=True,anchor="end")
    sc.text(x+15,y+h-118,"High-quality modules",6); sc.text(x+w-15,y+h-118,summ['HIGH_QUALITY_MODULE_N'],12,bold=True,anchor="end")
    mm=[r for r in rows("06_results/gateA/A2_MODULE_METRICS.tsv") if r['QUALITY']=='HIGH_QUALITY']; x,y,w,h=boxes[3]
    bar_values(sc,(x+28,y+26,w-38,h-55),[r['module'] for r in mm],[fnum(r['ACTIVATION_EFFECT']) for r in mm],[PAL['D1'] if r['PEAK_TIME']=='D1' else PAL['D7'] for r in mm],"Activation")
    ms=subset(rows("06_results/gateA/A2_MODULE_SCORES.tsv"),module="D7_M3"); x,y,w,h=boxes[4]; paired_lines(sc,(x+30,y+28,w-42,h-52),ms,'donor','timepoint','score',['Skin','D1','D7','D30'])
    x,y,w,h=boxes[5]; sc.rect(x+36,y+44,78,72,fill=PAL['D7'],stroke=None,rx=8); sc.text(x+75,y+85,"D7_M3",10,"#FFFFFF",True,"middle"); sc.text(x+75,y+67,"CORE100",8,"#FFFFFF",True,"middle"); sc.line(x+119,y+80,x+155,y+80,"#4A5560",1.2); sc.rect(x+157,y+55,48,50,fill="#EEF0F4",stroke="#6A3D9A",rx=5); sc.text(x+181,y+78,"LOCKED",7,PAL['keloid'],True,"middle"); sc.text(x+w/2,y+19,"100 locked genes • prespecified direction",6,anchor="middle")
    sc.export(MAIN/"Figure1")

def figure2():
    sc=Scene()
    boxes=[(18,335,150,275),(180,335,150,275),(342,335,150,275),(18,28,231,275),(261,28,231,275)]
    titles=["Discovery trajectory","Patient distributions","Early vs uninjured","Late vs early","Signature consistency"]
    for b,l,t in zip(boxes,"ABCDE",titles): sc.panel(*b,l,t)
    ms=subset(rows("06_results/gateA/A2_MODULE_SCORES.tsv"),module="D7_M3"); x,y,w,h=boxes[0]; paired_lines(sc,(x+30,y+34,w-42,h-62),ms,'donor','timepoint','score',['Skin','D1','D7','D30'])
    ps=rows("06_results/spectrum/GSE178411_PATIENT_STATE_TWPS.tsv"); x,y,w,h=boxes[1]; dot_groups(sc,(x+31,y+34,w-43,h-62),ps,'state','CORE100',['UNINJURED_SKIN','EARLY_WOUND','LATE_WOUND'])
    sr=rows("06_results/spectrum/GSE178411_TWPS_SUMMARY.tsv")
    for bi,label,note in [(2,'EARLY_VS_UNINJURED','Patient-aware P=1.11 × 10−20'),(3,'LATE_VS_EARLY','Patient-aware P=0.000993')]:
        r=subset(sr,contrast=label,signature='CORE100')[0]; x,y,w,h=boxes[bi]; forest(sc,(x+52,y+94,w-65,80),[label.replace('_',' ').title()],[fnum(r['difference'])],[fnum(r['ci_low'])],[fnum(r['ci_high'])],[PAL['D1'] if bi==2 else PAL['D30']],"Difference (95% CI)"); sc.wrapped(x+15,y+64,f"Difference {fnum(r['difference']):.3f}\n{note}",28,6,9)
    sig=[r for r in sr if r['contrast'] in ('EARLY_VS_UNINJURED','LATE_VS_EARLY')]; x,y,w,h=boxes[4]
    labels=['C25','C50','C100','Full']; order=['CORE25','CORE50','CORE100','FULL']; vals=[fnum(r['difference']) for r in sig]; lo=min(vals)-.15; hi=max(vals)+.15
    yy=axes(sc,x+28,y+45,w-42,h-78,lo,hi,labels,"Frozen difference")
    for i,s in enumerate(order):
        cx=x+28+(w-42)*(i+.5)/4
        e=subset(sig,contrast='EARLY_VS_UNINJURED',signature=s)[0]; l=subset(sig,contrast='LATE_VS_EARLY',signature=s)[0]
        sc.circle(cx-3,yy(fnum(e['difference'])),2.4,PAL['D1'],"#FFFFFF",.3); sc.circle(cx+3,yy(fnum(l['difference'])),2.4,PAL['D30'],"#FFFFFF",.3)
    sc.text(x+w/2,y+22,"Early−uninjured • Late−early",6.5,anchor="middle")
    sc.export(MAIN/"Figure2")

def figure3():
    sc=Scene()
    boxes=[(18,335,150,275),(180,335,150,275),(342,335,150,275),(18,28,150,275),(180,28,150,275),(342,28,150,275)]
    for b,l,t in zip(boxes,"ABCDEF",["Longitudinal design","Pre/post subjects","Delta comparison","Established keloid","Donor-level contrast","Locked sensitivities"]): sc.panel(*b,l,t)
    x,y,w,h=boxes[0]; sc.rect(x+20,y+155,47,45,fill=PAL['healthy'],stroke=None,rx=4); sc.text(x+43,y+178,"Healthy",7,"#FFFFFF",True,"middle"); sc.rect(x+83,y+155,52,45,fill=PAL['keloid-prone'],stroke=None,rx=4); sc.text(x+109,y+178,"Keloid-",6.2,"#FFFFFF",True,"middle"); sc.text(x+109,y+168,"prone",6.2,"#FFFFFF",True,"middle"); sc.line(x+43,y+147,x+43,y+100,"#66717B",1); sc.line(x+109,y+147,x+109,y+100,"#66717B",1); sc.text(x+w/2,y+76,"Baseline → post-wound",8,bold=True,anchor="middle"); sc.text(x+w/2,y+49,"n=4 healthy; n=8 keloid-prone",6,anchor="middle")
    gb=rows("06_results/gateB/corrected/GATE_B_TWPS_SCORES_CORRECTED.tsv"); x,y,w,h=boxes[1]
    vals=[fnum(r['CORE100']) for r in gb]; yy=axes(sc,x+30,y+35,w-43,h-65,min(vals)-.1,max(vals)+.1,['Pre','Post'],"CORE100")
    for ident in sorted(set(r['subject_id'] for r in gb)):
        rs={r['timepoint']:r for r in gb if r['subject_id']==ident}; pts=[]; col=PAL[rs['baseline']['group']]
        for i,t in enumerate(['baseline','post-wounding']): pts.append((x+30+(w-43)*(i+.5)/2,yy(fnum(rs[t]['CORE100']))))
        sc.polyline(pts,col,.8); [sc.circle(px,py,2,col,"#FFFFFF",.25) for px,py in pts]
    ch=rows("06_results/gateB/corrected/GATE_B_SUBJECT_CHANGES_CORRECTED.tsv"); x,y,w,h=boxes[2]; dot_groups(sc,(x+31,y+60,w-43,h-90),ch,'group','delta_TWPS',['healthy','keloid-prone'],ylabel='Δ TWPS'); sc.text(x+w/2,y+42,"Difference 0.094 • g=0.824",6.5,anchor="middle"); sc.text(x+w/2,y+29,"Exact permutation P=0.182",6.5,anchor="middle")
    gc=rows("06_results/gateC/GSE181316_DONOR_TWPS.tsv"); x,y,w,h=boxes[3]; dot_groups(sc,(x+31,y+48,w-43,h-78),gc,'group','CORE100',['normal_scar','keloid'],ylabel='CORE100'); sc.text(x+w/2,y+24,"3 donors per group",6,bold=True,anchor="middle")
    sm=rows("06_results/gateC/GSE181316_GATE_C_SUMMARY.tsv")[0]; x,y,w,h=boxes[4]; sc.rect(x+18,y+96,w-36,96,fill="#F5F1F8",stroke=None,rx=6); sc.text(x+w/2,y+166,"Hedges g = 2.310",9,bold=True,anchor="middle"); sc.text(x+w/2,y+145,"Cliff δ = 1 • superiority 9/9",6.7,anchor="middle"); sc.text(x+w/2,y+125,"Exact two-sided P = 0.10",7,bold=True,anchor="middle"); sc.text(x+w/2,y+57,"n=3 keloid vs 3 normal-scar donors",6,anchor="middle"); sc.text(x+w/2,y+39,"Precision limited; CI not emphasized",5.8,anchor="middle")
    gbs=rows("06_results/gateB/corrected/GATE_B_SIGNATURE_SENSITIVITY_CORRECTED.tsv"); gcs=rows("06_results/gateC/GSE181316_GATE_C_SUMMARY.tsv"); order=['CORE25','CORE50','CORE100','FULL']; x,y,w,h=boxes[5]
    vals=[fnum(subset(gbs,signature=s)[0]['difference']) for s in order]+[fnum(subset(gcs,signature=s)[0]['difference']) for s in order]
    yy=axes(sc,x+28,y+45,w-42,h-78,min(vals)-.12,max(vals)+.12,['C25','C50','C100','Full'],"Frozen difference")
    for i,s in enumerate(order):
        cx=x+28+(w-42)*(i+.5)/4; a=subset(gbs,signature=s)[0]; b=subset(gcs,signature=s)[0]
        sc.circle(cx-3,yy(fnum(a['difference'])),2.4,PAL['keloid-prone'],"#FFFFFF",.3); sc.circle(cx+3,yy(fnum(b['difference'])),2.4,PAL['keloid'],"#FFFFFF",.3)
    sc.text(x+w/2,y+22,"Keloid-prone • established keloid",5.9,bold=True,anchor="middle")
    sc.export(MAIN/"Figure3")

def figure4():
    sc=Scene()
    boxes=[(18,335,231,275),(261,335,231,275),(18,28,150,275),(180,28,150,275),(342,28,150,275)]
    for b,l,t in zip(boxes,"ABCDE",["HTS boundary distribution","HTS vs late","Spatial donor trajectories","Paired differences","Representative maps"]): sc.panel(*b,l,t)
    ps=rows("06_results/spectrum/GSE178411_PATIENT_STATE_TWPS.tsv"); x,y,w,h=boxes[0]; dot_groups(sc,(x+38,y+46,w-50,h-76),ps,'state','CORE100',['EARLY_WOUND','LATE_WOUND','HYPERTROPHIC_SCAR'])
    sr=subset(rows("06_results/spectrum/GSE178411_TWPS_SUMMARY.tsv"),contrast='HTS_VS_LATE',signature='CORE100')[0]; x,y,w,h=boxes[1]; forest(sc,(x+75,y+122,w-90,72),['HTS − late'],[fnum(sr['difference'])],[fnum(sr['ci_low'])],[fnum(sr['ci_high'])],[PAL['HYPERTROPHIC_SCAR']],"Difference (95% CI)"); sc.wrapped(x+18,y+91,"Difference −0.836; P=6.73 × 10−13",46,7,9); sc.rect(x+18,y+41,w-36,34,fill="#F7ECD4",stroke=None,rx=4); sc.text(x+w/2,y+61,"UNIVERSAL PATHOLOGICAL-SCAR",6.3,bold=True,anchor="middle"); sc.text(x+w/2,y+49,"PERSISTENCE NOT SUPPORTED",6.3,bold=True,anchor="middle")
    tr=subset(rows("06_results/spatial/GSE241124_DONOR_TRAJECTORIES.tsv"),signature='CORE100'); x,y,w,h=boxes[2]; paired_lines(sc,(x+30,y+48,w-42,h-78),tr,'donor','timepoint','TWPS_score',['SKIN','D1','D7','D30'],{'SKIN':PAL['Skin'],'D1':PAL['D1'],'D7':PAL['D7'],'D30':PAL['D30']}); sc.text(x+w/2,y+24,"4 donors; 2 overlap discovery",5.8,anchor="middle")
    pd=subset(rows("06_results/spatial/GSE241124_PAIRED_DONOR_DIFFERENCES.tsv"),signature='CORE100'); x,y,w,h=boxes[3]; vals1=[fnum(r['difference']) for r in pd if r['contrast']=='D7_VS_SKIN']; vals2=[fnum(r['difference']) for r in pd if r['contrast']=='D30_VS_D7']; allv=vals1+vals2; yy=axes(sc,x+30,y+58,w-42,h-88,min(allv)-.2,max(allv)+.2,['D7−Skin','D30−D7'],'Difference');
    for gi,valsx in enumerate([vals1,vals2]):
        cx=x+30+(w-42)*(gi+.5)/2
        for j,v in enumerate(valsx): sc.circle(cx+(j-1.5)*3,yy(v),2.2,PAL['D7'] if gi==0 else PAL['D30'],"#FFFFFF",.25)
    sc.text(x+w/2,y+37,"4/4 positive | 4/4 negative",5.8,bold=True,anchor="middle"); sc.text(x+w/2,y+24,"Exact P=0.125 for each",5.6,anchor="middle")
    spot=rows("06_results/spatial/GSE241124_SPOT_TWPS.tsv"); x,y,w,h=boxes[4]; sw=(w-24)/3
    for i,(s,t) in enumerate([('Donor1-Skin','Skin'),('Donor1-Wound7','D7'),('Donor1-Wound30','D30')]): spatial(sc,(x+8+i*sw,y+72,sw-5,125),spot,s,t)
    sc.text(x+w/2,y+48,"Donor1 selected by coordinate/image QC only",5.6,anchor="middle"); sc.text(x+w/2,y+34,"Per-slide display z; spots are descriptive",5.5,anchor="middle")
    sc.export(MAIN/"Figure4")

def figure5():
    sc=Scene(h=480)
    sc.panel(18,245,474,190,'A','Normal wound trajectory')
    stages=[('Skin',PAL['Skin'],0),('D1',PAL['D1'],1),('D7',PAL['D7'],1.2),('D30',PAL['D30'],.45)]
    xs=[72,190,308,426]; base=315
    for i,(lab,col,lev) in enumerate(stages):
        yy=base+lev*45; sc.circle(xs[i],yy,14,col,None,0); sc.text(xs[i],yy-2,lab,7,'#FFFFFF',True,'middle')
        if i<len(stages)-1: sc.line(xs[i]+16,yy,xs[i+1]-16,base+stages[i+1][2]*45,'#68727B',1.4)
    sc.text(255,274,"Activation followed by attenuation/remodeling; D30 is not assumed to be full resolution",6.5,anchor="middle")
    cards=[('Keloid-prone cohort','— Directionally larger post-injury TWPS change',PAL['keloid-prone']),('Established keloid','Higher observed TWPS than normal scar',PAL['keloid']),('HTS boundary','Universal persistence not supported',PAL['HYPERTROPHIC_SCAR']),('Genetics','No detectable program-level enrichment',PAL['Other'])]
    for i,(head,body,col) in enumerate(cards):
        x=18+i*119; sc.rect(x,115,109,92,fill=col,stroke=None,rx=6); sc.text(x+54.5,180,head,7,'#FFFFFF',True,'middle'); sc.wrapped(x+10,158,body,24,6.2,8,'#FFFFFF')
    sc.rect(55,36,400,48,fill="#EEF1F4",stroke=None,rx=8); sc.text(255,64,"TWPS captures a dynamic wound-response state with",8,bold=True,anchor="middle"); sc.text(255,50,"directional patterns in keloid-related states",8,bold=True,anchor="middle"); sc.text(255,20,"Association • attenuation • boundary testing — no causal or genetic-mediation claim",5.8,anchor="middle")
    sc.export(MAIN/"Figure5")

def supplementary():
    # S1 discovery robustness
    sc=Scene(title="Supplementary Figure S1 | GSE241132 discovery robustness")
    summ={r['metric']:fnum(r['value']) for r in rows("06_results/gateA/gate_a1_summary.tsv")}; a2={r['metric']:fnum(r['value']) for r in rows("06_results/gateA/gate_a2_summary.tsv")}
    sc.panel(18,335,231,275,'A','Transient tiers'); bar_values(sc,(48,380,180,185),['Tier1','Tier2','Tier3'],[summ['tier1_total'],summ['tier2_total'],summ['tier3_total']],[PAL['D1'],PAL['D7'],PAL['Other']],"Genes")
    sc.panel(261,335,231,275,'B','LODO stability'); ld=rows("06_results/gateA/LODO_STABILITY.tsv"); bar_values(sc,(300,390,170,170),[r['omitted_donor'][-2:]+' '+r['program'] for r in ld],[fnum(r['jaccard']) for r in ld],[PAL['D7'] if r['program']=='D7' else PAL['D1'] for r in ld],"Jaccard",False)
    sc.panel(18,28,474,275,'C','Fibroblast composition context'); fc=rows("06_results/gateA/FIBROBLAST_DONOR_TIME_COUNTS.tsv"); labels=[r['donor'][-2:]+'-'+r['timepoint'] for r in fc]; bar_values(sc,(50,72,420,185),labels,[fnum(r['n_fibroblast_cells']) for r in fc],[PAL.get(r['timepoint'],PAL['Other']) for r in fc],"Cells",True); sc.text(255,44,"Cell counts are technical context; biological n=3 donors",6,anchor="middle"); sc.export(SUPP/"FigureS1")
    # S2
    sc=Scene(title="Supplementary Figure S2 | Five-module robustness")
    mm=[r for r in rows("06_results/gateA/A2_MODULE_METRICS.tsv") if r['QUALITY']=='HIGH_QUALITY'];
    sc.panel(18,335,231,275,'A','Activation and attenuation'); labels=[r['module'] for r in mm]; bar_values(sc,(55,382,175,180),labels,[fnum(r['ACTIVATION_EFFECT']) for r in mm],[PAL['D1'] if r['PEAK_TIME']=='D1' else PAL['D7'] for r in mm],"Activation")
    sc.panel(261,335,231,275,'B','Attenuation'); bar_values(sc,(300,382,170,180),labels,[fnum(r['ATTENUATION_EFFECT']) for r in mm],[PAL['D1'] if r['PEAK_TIME']=='D1' else PAL['D7'] for r in mm],"Attenuation")
    st=rows("06_results/gateA/A2_SUBTYPE_MODULE_SENSITIVITY.tsv"); sc.panel(18,28,474,275,'C','Frozen fibroblast-subtype sensitivity');
    y=260
    for r in st: sc.text(42,y,f"{r['subtype']}  {r['module']}",6); sc.rect(170,y-3,220*fnum(r['activation'])/3,6,fill=PAL['D7'],stroke=None); sc.text(402,y,f"direction preserved: {r['direction_preserved_2of3']}",5.5); y-=18
    sc.export(SUPP/"FigureS2")
    # S3
    sc=Scene(title="Supplementary Figure S3 | GSE113619 locked-signature sensitivity")
    gb=rows("06_results/gateB/corrected/GATE_B_SIGNATURE_SENSITIVITY_CORRECTED.tsv"); order=['CORE25','CORE50','CORE100','FULL']; rr=[subset(gb,signature=s)[0] for s in order]
    sc.panel(18,335,474,275,'A','Corrected effect estimates'); forest(sc,(115,390,335,160),order,[fnum(r['hedges_g']) for r in rr],[fnum(r['ci_low']) for r in rr],[fnum(r['ci_high']) for r in rr],[PAL['keloid-prone']]*4,"Hedges g (95% CI)")
    sc.panel(18,28,474,275,'B','Exact permutation results'); bar_values(sc,(70,78,380,160),order,[fnum(r['permutation_p']) for r in rr],[PAL['Other']]*4,"P",False); sc.text(255,52,"CORE100 coverage 96/100; final corrected implementation",6.2,bold=True,anchor="middle"); sc.export(SUPP/"FigureS3")
    # S4
    sc=Scene(title="Supplementary Figure S4 | GSE181316 small-sample sensitivity")
    gc=rows("06_results/gateC/GSE181316_GATE_C_SUMMARY.tsv"); order=['CORE25','CORE50','CORE100','FULL']; rr=[subset(gc,signature=s)[0] for s in order]
    sc.panel(18,335,231,275,'A','All donor observations'); dot_groups(sc,(55,382,175,180),rows("06_results/gateC/GSE181316_DONOR_TWPS.tsv"),'group','CORE100',['normal_scar','keloid'])
    sc.panel(261,335,231,275,'B','Locked-signature direction'); bar_values(sc,(300,382,170,180),order,[fnum(r['difference']) for r in rr],[PAL['keloid']]*4,"Difference")
    sc.panel(18,28,474,275,'C','Precision warning'); sc.rect(65,95,380,130,fill="#F5F1F8",stroke=None,rx=8); sc.text(255,188,"All signatures: Cliff δ=1; pairwise superiority 9/9",8,bold=True,anchor="middle"); sc.text(255,160,"Exact two-sided P=0.10",8,anchor="middle"); sc.text(255,130,"n=3 per group — bootstrap standardized-effect intervals are unstable",7,bold=True,anchor="middle"); sc.export(SUPP/"FigureS4")
    # S5
    sc=Scene(title="Supplementary Figure S5 | Expanded GSE178411 validation")
    ps=rows("06_results/spectrum/GSE178411_PATIENT_STATE_TWPS.tsv"); groups=['UNINJURED_SKIN','EARLY_WOUND','LATE_WOUND','CHRONIC_WOUND','NORMAL_SCAR','HYPERTROPHIC_SCAR']
    sc.panel(18,335,474,275,'A','All frozen patient-state distributions'); dot_groups(sc,(55,385,415,175),ps,'state','CORE100',groups)
    sr=rows("06_results/spectrum/GSE178411_TWPS_SUMMARY.tsv"); rr=[r for r in sr if r['signature']=='CORE100']; sc.panel(18,28,231,275,'B','Primary contrasts'); forest(sc,(115,88,110,145),[r['contrast'].replace('_VS_',' vs ') for r in rr],[fnum(r['difference']) for r in rr],[fnum(r['ci_low']) for r in rr],[fnum(r['ci_high']) for r in rr],None,"Difference")
    sc.panel(261,28,231,275,'C','Continuous-time traceability'); sc.wrapped(280,220,"Frozen n=106, slope −0.151, 95% CI −0.183 to −0.120, P=1.63×10⁻¹⁴.",43,7,10); sc.wrapped(280,155,"The exact 106-row input was reconstructed; coefficient, CI and P matched the frozen result exactly.",43,7,10); sc.export(SUPP/"FigureS5")
    # S6 all spatial
    sc=Scene(h=690,title="Supplementary Figure S6 | All GSE241124 spatial maps")
    spot=rows("06_results/spatial/GSE241124_SPOT_TWPS.tsv"); samples=[]
    for d in ['Donor1','Donor2','Donor3','Donor4']:
        for t in [('Skin','Skin'),('Wound1','D1'),('Wound7','D7'),('Wound30','D30')]: samples.append((f'{d}-{t[0]}',f'{d} {t[1]}'))
    for i,(s,t) in enumerate(samples):
        col=i%4; row=i//4; x=30+col*119; y=505-row*145; spatial(sc,(x,y,92,100),spot,s,t)
    sc.text(255,20,"Per-slide display z-score; spots are not biological replicates",6.3,bold=True,anchor="middle"); sc.export(SUPP/"FigureS6")
    # S7 rho
    sc=Scene(title="Supplementary Figure S7 | Descriptive TWPS–fibroblast spatial association")
    ar=rows("06_results/spatial/GSE241124_TWPS_FIBROBLAST_ASSOCIATION.tsv"); sc.panel(18,150,474,460,'A','One Spearman rho per slide');
    times=['SKIN','D1','D7','D30']; yy=axes(sc,58,205,410,330,0,1,times,"Spearman rho")
    for i,t in enumerate(times):
        rs=[r for r in ar if r['timepoint']==t]; cx=58+410*(i+.5)/4
        for j,r in enumerate(rs): sc.circle(cx+(j-1.5)*7,yy(fnum(r['spearman_rho'])),3,PAL['Skin'] if t=='SKIN' else PAL.get(t,PAL['Other']),"#FFFFFF",.3)
    sc.text(255,120,"16/16 positive • median rho=0.506221",8,bold=True,anchor="middle"); sc.text(255,100,"Descriptive spatial association; no spot-level inference",6.3,anchor="middle"); sc.export(SUPP/"FigureS7")
    # S8 genetics
    sc=Scene(title="Supplementary Figure S8 | Prospective genetic boundary")
    pr=rows("06_results/genetics/CORE100_PRIMARY_MAGMA_RESULT.tsv"); sc.panel(18,335,231,275,'A','CORE100 competitive MAGMA'); forest(sc,(85,405,140,125),[r['ancestry'] for r in pr],[fnum(r['beta']) for r in pr],[fnum(r['beta'])-fnum(r['SE']) for r in pr],[fnum(r['beta'])+fnum(r['SE']) for r in pr],[PAL['Other']]*2,"Beta ± SE"); sc.text(133,370,"EUR P=0.902 • AFR P=0.783",6,anchor="middle")
    ms=rows("06_results/genetics/MAGMA_SIGNATURE_SENSITIVITY.tsv"); sc.panel(261,335,231,275,'B','Signature sensitivities');
    order=['CORE25','CORE50','CORE100','FULL']; x,y,w,h=(300,390,170,170); vals=[]; labs=[]
    for anc in ['EUR','AFR']:
        for sig in order:
            rr=subset(ms,ancestry=anc,signature=sig)
            if not rr and sig=='CORE100': rr=subset(pr,ancestry=anc)
            if rr: vals.append(fnum(rr[0]['beta'])); labs.append(anc[0]+sig.replace('CORE','').replace('FULL','F'))
    bar_values(sc,(x,y,w,h),labs,vals,[PAL['Other']]*len(vals),"Beta")
    sc.panel(18,28,474,275,'C','Secondary convergence summary');
    for i,(lab,val) in enumerate([('Recoverable FUMA subset overlap','0 / 115'),('Relevant-tissue GPGE CORE100','0'),('Colocalized CORE100','0')]): sc.text(65,225-i*55,lab,7); sc.text(430,225-i*55,val,10,bold=True,anchor="end")
    sc.rect(65,54,380,35,fill="#EEF1F4",stroke=None,rx=5); sc.text(255,74,"NO DETECTABLE PROGRAM-LEVEL GENETIC SUPPORT",7,bold=True,anchor="middle"); sc.export(SUPP/"FigureS8")

def copy_tables_and_source():
    mapping={
      'TableS1_Dataset_Roles_Biological_N.tsv':'06_results/final_audit/FINAL_DATASET_LEDGER.tsv',
      'TableS2_CORE100_Genes.tsv':'06_results/gateA/TWPS_PRIMARY_D7_M3_CORE100.tsv',
      'TableS3_Module_Gene_Membership.tsv':'06_results/gateA/A2_MODULE_MEMBERSHIP.tsv',
      'TableS4_GateB_Corrected_Statistics.tsv':'06_results/gateB/corrected/GATE_B_SIGNATURE_SENSITIVITY_CORRECTED.tsv',
      'TableS5_GateC_Statistics.tsv':'06_results/gateC/GSE181316_GATE_C_SUMMARY.tsv',
      'TableS6_GSE178411_Validation_Statistics.tsv':'06_results/spectrum/GSE178411_TWPS_SUMMARY.tsv',
      'TableS7_GSE241124_Spatial_Statistics.tsv':'06_results/spatial/GSE241124_SPATIAL_TWPS_SUMMARY.tsv',
      'TableS8_MAGMA_Genetic_Results.tsv':'06_results/genetics/MAGMA_SIGNATURE_SENSITIVITY.tsv',
      'TableS9_Deferred_Excluded_Dataset_Provenance.tsv':'06_results/final_audit/NONFINAL_BRANCH_AUDIT.tsv',
      'TableS10_Software_Reproducibility.tsv':'06_results/final_audit/REPRODUCIBILITY_MANIFEST.tsv'}
    for out,src in mapping.items(): shutil.copyfile(ROOT/src,TABLES/out)
    # Long-format main source-data packages, copied from frozen rows without recomputation.
    specs={
      'Figure1_SourceData.tsv':[('B','GSE241132','06_results/gateA/PSEUDOBULK_SAMPLE_MAP.tsv'),('C_D','GSE241132','06_results/gateA/A2_MODULE_METRICS.tsv'),('E','GSE241132','06_results/gateA/A2_MODULE_SCORES.tsv'),('F','GSE241132','06_results/gateA/TWPS_PRIMARY_D7_M3_CORE100.tsv')],
      'Figure2_SourceData.tsv':[('A','GSE241132','06_results/gateA/A2_MODULE_SCORES.tsv'),('B','GSE178411','06_results/spectrum/GSE178411_PATIENT_STATE_TWPS.tsv'),('C_D_F','GSE178411','06_results/spectrum/GSE178411_TWPS_SUMMARY.tsv')],
      'Figure3_SourceData.tsv':[('B','GSE113619','06_results/gateB/corrected/GATE_B_TWPS_SCORES_CORRECTED.tsv'),('C','GSE113619','06_results/gateB/corrected/GATE_B_SUBJECT_CHANGES_CORRECTED.tsv'),('D_E','GSE181316','06_results/gateC/GSE181316_DONOR_TWPS.tsv'),('F','validation','06_results/gateC/GSE181316_SIGNATURE_SENSITIVITY.tsv')],
      'Figure4_SourceData.tsv':[('A','GSE178411','06_results/spectrum/GSE178411_PATIENT_STATE_TWPS.tsv'),('B','GSE178411','06_results/spectrum/GSE178411_TWPS_SUMMARY.tsv'),('C_D_E','GSE241124','06_results/spatial/GSE241124_PAIRED_DONOR_DIFFERENCES.tsv'),('F','GSE241124','06_results/spatial/GSE241124_SPOT_TWPS.tsv')],
      'Figure5_SourceData.tsv':[('A','synthesis','06_results/final_audit/CLAIM_STRENGTH_MATRIX.tsv')]}
    for out,groups in specs.items():
        records=[]; allcols=set()
        for panel,dataset,src in groups:
            rr=rows(src)
            for r in rr:
                z={'panel':panel,'dataset':dataset,'source_file':src}; z.update(r); records.append(z); allcols.update(z)
        cols=['panel','dataset','source_file']+sorted(allcols-{'panel','dataset','source_file'})
        with open(SOURCE/out,'w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',extrasaction='ignore'); w.writeheader(); w.writerows(records)
    with open(SOURCE/'Supplementary_SourceData_Manifest.tsv','w',newline='',encoding='utf-8') as f:
        w=csv.writer(f,delimiter='\t'); w.writerow(['figure','source_files'])
        w.writerows([
          ['FigureS1','gate_a1_summary.tsv;LODO_STABILITY.tsv;FIBROBLAST_DONOR_TIME_COUNTS.tsv'],
          ['FigureS2','A2_MODULE_METRICS.tsv;A2_SUBTYPE_MODULE_SENSITIVITY.tsv'],
          ['FigureS3','GATE_B_SIGNATURE_SENSITIVITY_CORRECTED.tsv'],['FigureS4','GSE181316_GATE_C_SUMMARY.tsv;GSE181316_DONOR_TWPS.tsv'],
          ['FigureS5','GSE178411_PATIENT_STATE_TWPS.tsv;GSE178411_TWPS_SUMMARY.tsv'],['FigureS6','GSE241124_SPOT_TWPS.tsv'],
          ['FigureS7','GSE241124_TWPS_FIBROBLAST_ASSOCIATION.tsv'],['FigureS8','CORE100_PRIMARY_MAGMA_RESULT.tsv;MAGMA_SIGNATURE_SENSITIVITY.tsv']])

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--main-only',action='store_true',help='Render only main Figures 1–5 from frozen source tables.')
    args=parser.parse_args()
    for d in [MAIN,SUPP,TABLES,SOURCE]: d.mkdir(parents=True,exist_ok=True)
    figure1(); figure2(); figure3(); figure4(); figure5()
    if not args.main_only:
        supplementary(); copy_tables_and_source()
    print('PHASE10B_RENDER_COMPLETE')

if __name__=='__main__': main()
