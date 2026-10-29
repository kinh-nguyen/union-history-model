import sys, os, pandas as pd, numpy as np
# Build the grouped post-marriage data for one or more countries: prep_post.py CC [CC ...]
# CC8 names the cell container, so the same script serves the 37 sub-Saharan countries
# out of cc8.csv.bz2 and the 33 extension countries out of cc8_global.csv.bz2.
# Grouping key is (age at marriage, duration, cell, cohort, survey) -- deliberately NOT age at
# debut, which the post-marriage rates do not depend on; keeping it would cost 3.75x the chains.
so = ["V","X","M","D1","W1","R","D2","W2"]
cc = pd.read_csv(os.environ.get("CC8", "../MultistageSurv/data/cc8.csv.bz2"))
for c in (sys.argv[1:] or ["RW", "ST"]):
    d = cc[(cc.cc == c) & (cc.sex == "f") & (cc.A == "M")].copy()
    assert len(d), f"no female M-episodes for {c}"
    d["start"] = d.start.where(d.start != 0, 1); d["dur"] = d.end - d.start; d = d[d.dur >= 1]
    d["Zi"] = d.Z.map({s: i for i, s in enumerate(so)}) - 2; d["coh"] = (d.yob - 1975)/10
    g = d.groupby(["start","dur","Zi","coh","sv"], as_index=False).n.sum()
    np.savez(f"post_{c}.npz", start=g.start.values, dur=g.dur.values, Z=g.Zi.values, coh=g.coh.values, n=g.n.values, sv=g.sv.astype(str).values)
    print(c, "rows", len(g), "weight", round(g.n.sum()), "surveys", g.sv.nunique(), g.groupby("Zi").n.sum().round().to_dict())
