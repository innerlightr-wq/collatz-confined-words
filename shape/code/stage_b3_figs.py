import sys, csv, math, collections, statistics
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
def load(p):
    out=[]
    for d in csv.DictReader(open(p)):
        d["r"]=int(d["r"]); d["theta"]=float(d["theta"]); d["D"]=int(d["D"])
        d["margin"]=float(d["margin"]); d["tied"]=d["tied"]=="True"
        d["S4"]=float(d["S4"]) if d.get("S4") else None
        d["delta"]=int(d["delta"]) if d.get("delta") else None
        out.append(d)
    return out
E=load("data/balance_explore.csv"); H=load("data/balance_holdout.csv")
fig,ax=plt.subplots(2,2,figsize=(12.5,9.5))
C0,C1="#2f6f9f","#c4432b"
# (a) margin vs r
a=ax[0,0]
fin=[d for d in E if math.isfinite(d["margin"]) and d["margin"]>0]
a.scatter([d["r"] for d in fin],[d["margin"] for d in fin],s=1,alpha=.12,c=C0,label="untied (m>0)",rasterized=True)
t=[d for d in E if d["tied"]]
a.scatter([d["r"] for d in t],[1e-7]*len(t),s=6,c=C1,label="tie (m = 0 exactly)")
RB=[1,5,10,15,20,30,50,100,400]
mx=[];my=[]
for i in range(len(RB)-1):
    s_=[d["margin"] for d in E if RB[i]<=d["r"]<RB[i+1] and math.isfinite(d["margin"])]
    if s_: mx.append((RB[i]+RB[i+1])/2); my.append(statistics.median(s_))
a.plot(mx,my,c="k",lw=2,label="median m")
a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel("r"); a.set_ylabel("balance margin m")
a.set_title("(a) the profile gets MORE balanced as r grows\nyet ties stop happening"); a.legend(fontsize=8)
# (b) median m vs median |Delta|
a=ax[0,1]
bx,bm,bd,bt=[],[],[],[]
for i in range(len(RB)-1):
    s_=[d for d in H if RB[i]<=d["r"]<RB[i+1] and d["delta"] is not None and math.isfinite(d["margin"])]
    if len(s_)<20: continue
    bx.append((RB[i]+RB[i+1])/2); bm.append(statistics.median([d["margin"] for d in s_]))
    bd.append(max(1,statistics.median([abs(d["delta"]) for d in s_]))); bt.append(sum(d["tied"] for d in s_)/len(s_))
a.plot(bx,bm,"o-",c=C0,label="median relative margin m")
a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel("r"); a.set_ylabel("m",color=C0)
a2=a.twinx(); a2.plot(bx,bd,"s-",c=C1,label="median |Delta| (exact integer gap)")
a2.set_yscale("log"); a2.set_ylabel("|Delta|",color=C1)
a.set_title("(b) HOLDOUT: relative balance improves, absolute\ninteger gap explodes -- ties follow |Delta|")
a.legend(fontsize=8,loc="center left"); a2.legend(fontsize=8,loc="lower right")
# (c) S4 vs theta at matched r
a=ax[1,0]
for lo,hi,c,mk in ((10,20,"#2f6f9f","o"),(20,40,"#4a9a5a","s"),(40,80,"#c4432b","^")):
    pts=[(d["theta"],d["S4"]) for d in E if lo<=d["r"]<hi and d["S4"] is not None]
    if not pts: continue
    TH=np.arange(0.13,0.91,0.08); xs=[];ys=[]
    for i in range(len(TH)-1):
        s_=[p[1] for p in pts if TH[i]<=p[0]<TH[i+1]]
        if s_: xs.append((TH[i]+TH[i+1])/2); ys.append(statistics.median(s_))
    a.plot(xs,ys,mk+"-",c=c,label="r %d-%d"%(lo,hi))
a.set_xlabel("theta"); a.set_ylabel("median S4 (TV to best ideal kernel)")
a.set_title("(c) distance to the ideal kernel grows with theta\n(P3's direction: holds)"); a.legend(fontsize=8)
# (d) m vs D within slopes
a=ax[1,1]
per=collections.defaultdict(list)
for d in E:
    if math.isfinite(d["margin"]) and d["r"]>=20: per[d["slope"]].append(d)
cs=[]
for n,v in per.items():
    if len(v)<50: continue
    x=np.array([d["D"] for d in v]); y=np.array([d["margin"] for d in v])
    rx=np.argsort(np.argsort(x)); ry=np.argsort(np.argsort(y))
    cs.append(float(np.corrcoef(rx,ry)[0,1]))
a.hist(cs,bins=20,color=C0,edgecolor="k",lw=.5)
a.axvline(0,c="k",lw=1.5); a.axvline(statistics.median(cs),c=C1,lw=2,ls="--",
                                     label="median %+.3f"%statistics.median(cs))
a.set_xlabel("within-slope Spearman(m, D)  [r >= 20]"); a.set_ylabel("slopes")
a.set_title("(d) the mechanism needs this clearly positive;\nit is centred on zero (EXPLORE, r >= 20)"); a.legend(fontsize=8)
plt.tight_layout(); plt.savefig("figs/balance_summary.png",dpi=130); print("wrote figs/balance_summary.png")
# holdout figure
fig,ax=plt.subplots(1,2,figsize=(11,4.2))
a=ax[0]
fin=[d for d in H if math.isfinite(d["margin"]) and d["margin"]>0]
a.scatter([d["r"] for d in fin],[d["margin"] for d in fin],s=2,alpha=.2,c=C0,rasterized=True)
t=[d for d in H if d["tied"]]
a.scatter([d["r"] for d in t],[1e-7]*len(t),s=14,c=C1,label="tie (n=%d)"%len(t))
a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel("r"); a.set_ylabel("m")
a.set_title("HOLDOUT: ties only at r < 10"); a.legend(fontsize=8)
a=ax[1]
bx,bt2=[],[]
RB2=[1,5,10,15,20,30,50,100,400]
for i in range(len(RB2)-1):
    s_=[d for d in H if RB2[i]<=d["r"]<RB2[i+1] and math.isfinite(d["margin"])]
    if len(s_)<20: continue
    s2=sorted(s_,key=lambda d:d["margin"]); top=s2[:max(1,len(s2)//10)]
    bx.append((RB2[i]+RB2[i+1])/2); bt2.append(sum(d["tied"] for d in top)/len(top))
a.plot(bx,bt2,"o-",c=C1); a.set_xscale("log"); a.set_xlabel("r")
a.set_ylabel("tie rate among the 10% most balanced rows")
a.set_title("HOLDOUT: maximal balance, zero ties, once r >= 10")
plt.tight_layout(); plt.savefig("figs/balance_holdout.png",dpi=130); print("wrote figs/balance_holdout.png")
