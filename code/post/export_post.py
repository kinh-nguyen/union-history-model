# Write the grouped post-marriage data and the two bases as CSV for the TMB driver,
# using exactly the grouping and bases of the JAX prototype (post.py).
import sys, numpy as np, pandas as pd
sys.path.insert(0, ".."); import post as PM
for cc in sys.argv[1:]:
    D = PM.load(cc, "..")
    st, coh = np.asarray(D["start"]), np.asarray(D["coh"])
    rows = pd.DataFrame(dict(g=D["g"], dur=D["dur"], Z=D["Z"], n=np.asarray(D["n"]))).sort_values(["g", "dur"])
    G = len(st); dmax = rows.groupby("g").dur.max().reindex(range(G)).fillna(0).astype(int).values
    roff = np.r_[0, np.cumsum(rows.groupby("g").size().reindex(range(G)).fillna(0).astype(int).values)]
    gy0 = np.round(1975 + 10*coh + st).astype(int) - PM.Y0
    pd.DataFrame(dict(start=st, coh=coh, dmax=dmax, gy0=gy0, ga0=np.round(st).astype(int))).to_csv(f"groups_{cc}.csv", index=False)
    pd.DataFrame(dict(off=roff)).to_csv(f"roff_{cc}.csv", index=False)
    rows.to_csv(f"rows_{cc}.csv", index=False)
    print(cc, G, len(rows), dmax.max())
np.savetxt("PB.csv", np.asarray(PM.PB), delimiter=","); np.savetxt("AB.csv", np.asarray(PM.AB), delimiter=",")
