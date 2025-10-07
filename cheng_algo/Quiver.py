import numpy as np
from time import time,sleep
import matplotlib.pyplot as plt
import networkx as nx
import Field 

## QUIVER REP
#definition
class Quiver:
    def __init__(self, vertices, edges,Ve,field=Field.Field('R'),grid=False):
        self.vertices=vertices#Q0
        self.edges = edges 
        self.Ve=Ve
        self.field=field
        self.grid=grid
        self.Ve_rounded=list(self.Ve.values())
        if self.field.descr in ['R','C']:
            self.Ve_rounded=[np.round(self.Ve[ae],2) for ae in self.Ve.keys()]
        if self.field.descr in ['Q']:
            self.Ve_rounded=[  np.array([np.round( float(f),2) for f in self.Ve[ae].flatten() ]).reshape(self.Ve[ae].shape)  for ae in self.Ve.keys()]
        #sort edges by starting/arriving extremity
        edges_in=dict(zip(vertices,[[] for _ in vertices]))
        edges_out=dict(zip(vertices,[[] for _ in vertices]))
        for e in edges.keys():
            edges_in[edges[e][1]].append(e)
            edges_out[edges[e][0]].append(e)
        self.edges_in,self.edges_out= edges_in,edges_out
        self.spaces={}
        for x in vertices:
            dims_x= [self.Ve[e].shape[0] for e in self.edges_in[x]]+[self.Ve[e].shape[1] for e in self.edges_out[x]]
            try:
                assert(len(set(dims_x))<=1)
            except:
                raise ValueError("dims of maps do not coincinde at vertex "+str(x))
                print(x,dims_x)
                print( [(e,self.Ve[e]) for e in self.edges_in[x]])
                print([(e,self.Ve[e]) for e in self.edges_out[x]])
            self.spaces[x]= dims_x[0]
            
    def is_zero(self):
        return not(any(self.spaces.values()))
    
    def display_graph(self,label="",verbose=True):
        plt.ion()
        if verbose:
            ax = plt.gca()
            ax.set_title(label)
            G = nx.DiGraph()
            G.add_nodes_from(self.vertices)
            edges= [(e[0],e[1]) for e in self.edges.values()]
            G.add_edges_from(edges)
            pos=nx.spring_layout(G,seed=7)
            if self.grid:
                pos = dict(zip(self.vertices,[np.asarray(v) for v in self.vertices ]))
            nx.draw(G,pos,ax=ax)
            edge_lab=dict(zip(edges, self.Ve_rounded))        
            nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_lab, font_color='red')
            nx.draw_networkx_labels(G, pos,labels=self.spaces, font_color='red')
            _ = ax.axis('off')
        plt.draw()
        plt.show()
        plt.pause(0.05)
  #display
    def __str__(self):
        return "vertices="+str(self.vertices)+" edges="+str(self.edges)+" Ve="+str(self.Ve_rounded)   
      
    def __repr__(self):
        return str(self)
    
    