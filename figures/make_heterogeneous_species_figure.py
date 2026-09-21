import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).parents[1]
geo=json.loads((ROOT/'..'/'interface_liver_twin'/'true_3d_voronoi.json').read_text())['cells']
rows=json.loads((ROOT/'..'/'interface_liver_twin'/'direct_sbml_spatial_results.json').read_text())['rows']
by={}
for r in rows: by.setdefault(int(r['path_id'].split('_')[0][1:]),[]).append(r)
pts=[]; vals={k:[] for k in ('oxygen','apap','napqi')}
for c in geo:
    rr=by.get(c['id'],[])
    if not rr: continue
    pts.append(c['center'])
    for k in vals: vals[k].append(float(np.mean([x[k] for x in rr])))
p=np.asarray(pts)
fig=plt.figure(figsize=(13,4.5),dpi=220)
for i,(k,label,cmap) in enumerate([('oxygen','Oxygen','viridis'),('apap','APAP','Blues'),('napqi','NAPQI','magma')],1):
    ax=fig.add_subplot(1,3,i,projection='3d'); v=np.asarray(vals[k]); sc=ax.scatter(p[:,0],p[:,1],p[:,2],c=v,s=34,cmap=cmap,alpha=.9,depthshade=False)
    ax.set_title(label); ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([]); ax.set_box_aspect((1,1,1)); ax.view_init(24,35); fig.colorbar(sc,ax=ax,shrink=.65,pad=.02)
fig.suptitle('Heterogeneous physiological fields on the virtual Voronoi tissue')
fig.tight_layout(); fig.savefig(ROOT/'heterogeneous_species_3d.png',facecolor='white'); fig.savefig(ROOT/'heterogeneous_species_3d.svg',facecolor='white')
print('saved',len(p),'lobule fields')
