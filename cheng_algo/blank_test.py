<<<<<<< HEAD
#%%
import numpy as np
import auxHN as aux
import copy
import Field as field
import matplotlib.pyplot as plt
from scipy import stats
from time import sleep
import mainHN as main
#aezaeaz

'''
V=aux.random_grid_indec([3,4], 5,3,3)

V.display_graph("V")


HNV=main.computeHN_sub(V,(0,0),verbose=False)

print(HNV)
plt.show()



'''
field=field.Field('F_2')



V1,supp1= aux.int_module_in_grid_quiver( [3,3], [(0,0)],[(2,0),(0,2),(1,1)],field)
V2,supp2= aux.int_module_in_grid_quiver( [3,3], [(0,0)],[(1,0),(0,1)],field)
=======
import numpy as np
import auxHN as aux
import copy
import matplotlib.pyplot as plt
from scipy import stats
#import mainHN as main
#aezaeaz


V=aux.random_grid_indec([5,5], 4)

#V.display_graph("V")

'''

V1,supp1= aux.int_module_in_grid_quiver( [2,2], [(0,0)],[(1,1)],field)
V2,supp2= aux.int_module_in_grid_quiver( [2,2], [(0,0)],[(1,0),(0,1)],field)
>>>>>>> main

x=(0,0)

V3=aux.ind_vertex(V1.vertices,V1.edges,x,field)
 
W= aux.direct_sum(aux.direct_sum(V1,V2),V3)


W.display_graph("1_x")


HNW=main.computeHN_sub(W,(0,0),verbose=True)

<<<<<<< HEAD
print("test:",HNW)

plt.show()

# %%
=======
print(HNW)
'''
>>>>>>> main
