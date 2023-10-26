
import itertools

import h5py as h5

import numpy as np
import plotly.graph_objects as go

with h5.File("./logentropy/output.h5", 'r') as f:

    n = f['n'][()]
    m = f['m'][()]
    l = f['l'][()]

    X = f['Yq'][()].T
    Y = f['logtemp'][()].T
    Z = f['lognb'][()].T

    f_in = f['f_in'][()].T
    f_out = f['f_out'][()].T

layout = go.Layout(width = 700, height =700)

fig = go.Figure(data=[go.Surface(x = Y[0,:,:], y = Z[0,:,:], z=f_out[0,:,:], colorscale = 'Blues')], layout=layout)
 
fig.add_scatter3d(x=Y[0,:,:].flatten(), y=Z[0,:,:].flatten(), z = f_in[0,:,:].flatten(), mode='markers', marker=dict(size=5, colorscale='Reds'))
 
fig.show()