import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection
ROOT=Path(__file__).parents[1]
geo=json.loads((ROOT/'..'/'interface_liver_twin'/'true_3d_voronoi.json').read_text())
net=json.loads((ROOT/'..'/'interface_liver_twin'/'lobule_sinusoid_network.json').read_text())
fig=plt.figure(figsize=(8,7),dpi=220); ax=fig.add_subplot(111,projection='3d')
segments=[]
for c in geo['cells']:
    if max(abs(x) for x in c['center'])>1.9: continue
    for face in c['faces']:
        for i in range(len(face)):
            a=c['vertices'][face[i]]; b=c['vertices'][face[(i+1)%len(face)]]
            segments.append([a,b])
ax.add_collection3d(Line3DCollection(segments,colors='#8aa0ad',linewidths=.35,alpha=.35))
cv=np.array([x['position'] for x in net['central_veins'] if max(abs(v) for v in x['position'])<1.9])
pt=np.array([x['position'] for x in net['portal_sites'] if max(abs(v) for v in x['position'])<1.9])
ax.scatter(cv[:,0],cv[:,1],cv[:,2],s=13,c='#c62828',label='central-vein generators',depthshade=False)
ax.scatter(pt[:,0],pt[:,1],pt[:,2],s=5,c='#c58a22',label='candidate portal sites',depthshade=False,alpha=.75)
ax.set(xlabel='x',ylabel='y',zlabel='z',title='Synthetic 3-D Voronoi virtual hepatic tissue')
ax.legend(loc='upper left',fontsize=8,frameon=False); ax.view_init(24,35); ax.set_box_aspect((1,1,1)); fig.tight_layout()
fig.savefig(ROOT/'voronoi_geometry_3d.png',facecolor='white'); fig.savefig(ROOT/'voronoi_geometry_3d.svg',facecolor='white')
print('saved geometry figure',len(segments),'edges',len(cv),'CVs',len(pt),'portal candidates')
