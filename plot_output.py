
import itertools

import numpy as np
import plotly.graph_objects as go

with open("output.dat", 'r') as f:
    lines = f.read().split("\n")

n, m, l = [int(i) for i in list(itertools.filterfalse(lambda x: x == '', lines[0].split(" ")))]
X = np.array( [float(i) for i in list(itertools.filterfalse(lambda x: x == '', lines[1].split(" ")))] )
Y = np.array( [float(i) for i in list(itertools.filterfalse(lambda x: x == '', lines[2].split(" ")))] )
Z = np.array( [float(i) for i in list(itertools.filterfalse(lambda x: x == '', lines[3].split(" ")))] )
f_in = np.array( [float(i) for i in list(itertools.filterfalse(lambda x: x == '', lines[4].split(" ")))] )
f_out = np.array( [float(i) for i in list(itertools.filterfalse(lambda x: x == '', lines[5].split(" ")))] )

X = X.reshape(n, m, l)
Y = Y.reshape(n, m, l)
Z = Z.reshape(n, m, l)
f_in = f_in.reshape(n, m, l)
f_out = f_out.reshape(n, m, l)


layout = go.Layout(width = 700, height =700)
 
fig = go.Figure(data=[go.Surface(x = X[:,:,0], y = Y[:,:,0], z=f_in[:,:,0], colorscale = 'Blues')], layout=layout)
 
fig.add_scatter3d(x=X[:,:,0].flatten(), y=Y[:,:,0].flatten(), z = f_out[:,:,0].flatten(), mode='markers', marker=dict(size=5, colorscale='Reds'))
 
fig.show()