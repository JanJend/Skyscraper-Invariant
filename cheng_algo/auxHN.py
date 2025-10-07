"""
Created on Thu Aug 12 19:50:59 2021

@author: marcf
"""


import numpy as np
import scipy.linalg
from time import time,sleep
import matplotlib.pyplot as plt
import traceback as tr
import sys,copy
from cfractions import Fraction

import Field 
import Quiver





#maximum computation error
epsilon=1e-10

#interval module of given support  in a  grid of dim xmax 
def int_module(vertices,edges,support,field):
    Ve={}
    for e in edges:
        Ve[e]= field.to_Field(np.ones((int(support[edges[e][1]]) ,int(support[edges[e][0]]))))
    return Quiver.Quiver(vertices,edges,Ve,field,grid=True)
    
def int_module_in_grid_quiver(xmax,sources,targets,field=Field.Field ("Q")):
    #build vertices and edges
    vertices=[(ij % xmax[0],ij//xmax[0]) for ij in range(xmax[0]*xmax[1])]
    edges={}
    for i in range(xmax[0]):
        for j in range(xmax[1]):
            if j<xmax[1]-1:
                edges[(i,j,1)]=[(i,j),(i,j+1)]
            if i<xmax[0]-1:
                edges[(i,j,0)]=[(i,j),(i+1,j)]
    #build support
    support={}
    for v in vertices:
        support[v]= any([all(np.array(v)>= np.array(s)) for s in sources]) and any([all(np.array(v)<= np.array(t)) for t in targets])
    return int_module(vertices,edges,support,field),support

def grid_indec(xmax, field=Field.Field("Q")):
    vertices=[(ij % xmax[0],ij//xmax[0]) for ij in range(xmax[0]*xmax[1])]
    edges={}
    Ve={}
    for i in range(xmax[0]):
        for j in range(xmax[1]):
            if j<xmax[1]-1:
                edges[(i,j,1)]=[(i,j),(i,j+1)]                   
            if i<xmax[0]-1:
                edges[(i,j,0)]=[(i,j),(i+1,j)]
    for e in edges.keys():
       if tuple(np.array(e)[:2])== (0,0):
           Ve[e]= field.to_Field(np.eye(2))
       elif e in [(1,0,1), (0,1,0)]:
           Ve[e]= field.to_Field(np.ones((1, 2)))
       elif e == (0,1,1):
           Ve[e]= np.array([[field.zero,field.one]])
       elif e == (1,0,0):
           Ve[e]= np.array([[field.one,field.zero]])    
       elif tuple(np.array(e)[:2]) in [ (2,0),(1,1), (0,2)]:
           Ve[e]= field.to_Field(np.zeros((0,1),dtype="i"))
       else:
           Ve[e]= field.to_Field(np.zeros((0, 0),dtype="i"))
    a=0
    return Quiver.Quiver(vertices,edges,Ve,field,grid=True)
    
def random_grid_indec(xmax, d_x,n_maps, d_y, field=Field.Field("Q")):
    vertices=[(ij % xmax[0],ij//xmax[0]) for ij in range(xmax[0]*xmax[1])]
    edges={}
    Ve={}
    for i in range(xmax[0]):
        for j in range(xmax[1]):
            if j<xmax[1]-1:
                edges[(i,j,1)]=[(i,j),(i,j+1)]                   
            if i<xmax[0]-1:
                edges[(i,j,0)]=[(i,j),(i+1,j)]
    # random proj on hyperplanes
    A=[np.random.randint(-5,5,size=(np.random.randint(1,d_y+1),d_x)) for _ in range(n_maps)]
    for e in edges.keys():
        if e[0]+e[1]<n_maps-2:
            Ve[e]= field.to_Field(np.eye(d_x))
        elif e[0]+e[1]==n_maps-1:
            Ve[e]= field.to_Field(np.zeros((0, np.shape(A[e[1]])[0]),dtype="i"))
        elif e[0]+e[1]==n_maps-2:
            Ve[e]= field.to_Field(A[e[1]+e[2]])
        else:
           Ve[e]= field.to_Field(np.zeros((0, 0),dtype="i"))
    return Quiver.Quiver(vertices,edges,Ve,field,grid=True)
    
def star_quiver(d_x,n_maps,d_y,field=Field.Field("Q")):
    vertices = list(range(n_maps+1))
    edges={}
    Ve={}
    A=[None]+[np.random.randint(-5,5,size=(np.random.randint(1,d_y+1),d_x)) for _ in range(n_maps)]
    for i in range(1,n_maps+1):
        edges[(0,i)]=[0,i]
        Ve[(0,i)]= field.to_Field(A[i])
    return Quiver.Quiver(vertices,edges,Ve,field,grid=False)
    
    
def ind_vertex(vertices,edges, x,field):
    Ve={}
    for i_e,e in edges.items():
        if e[0]== x:
            Ve[i_e]= field.to_Field(np.zeros((0,1),dtype="i"))
        elif e[1]==x:
            Ve[i_e]= field.to_Field(np.zeros((1,0),dtype="i"))
        else:
            Ve[i_e]=field.to_Field(np.zeros((0,0),dtype="i"))
    return Quiver.Quiver(vertices,edges,Ve,field)

def direct_sum(V,W):
    assert (V.edges.items() == W.edges.items() and V.field.descr==W.field.descr )
    maps=dict(zip(V.edges,[Field.build_block_diag(V.Ve[e],W.Ve[e],V.field) for e in V.edges]))
    Q=Quiver.Quiver(V.vertices, V.edges,maps,V.field,grid=V.grid)
    return Q
        
def subrep(V,Wx):
    assert (V.vertices==list(Wx.keys()))
    for x in V.vertices:
        try:
            Wx[x]=extract_basis(Wx[x], V.field)
        except:
            print("e",Wx[x])
    We=dict.copy(V.Ve)
    for e in V.edges.keys():
        VeWse=np.dot(We[e],Wx[V.edges[e][0]])
        We[e]=inverse_image_vect(Wx[V.edges[e][1]],VeWse, V.field)  
        '''
        try:
            We[e]=inverse_image_vect(Wx[V.edges[e][1]],VeWse, V.field)  
        except:
            print(repr(Wx[V.edges[e][0]])+" "+repr(Wx[V.edges[e][1]])+" "+repr(We[e]))
            raise ValueError("e1")
        '''
    return Quiver.Quiver(V.vertices, V.edges, We,V.field,grid=V.grid)
    
def quotientrep(V,Wx):
    assert (V.vertices==list(Wx.keys()))
    field=V.field
    #orthonormal bases
    WxT=dict.copy(Wx)
    for x in V.vertices:
        WxT[x]=complete_basis(Wx[x], field)
    #maps
    We=dict.copy(V.Ve)
    for e in V.edges.keys():
        VeWseT=np.dot(We[e],WxT[V.edges[e][0]])
        Vte_adpt=np.concatenate((WxT[V.edges[e][1]],Wx[V.edges[e][1]]),axis=1)
        We[e]=inverse_image_vect(Vte_adpt,VeWseT, V.field)[:WxT[V.edges[e][1]].shape[1]]   
    return Quiver.Quiver(V.vertices, V.edges, We,V.field,grid=V.grid)
   
 
## LINEAR ALGEBRA
    

# Row echelon form (Gaussion pivot)
<<<<<<< HEAD
def row_echelon(M_input,field, max_col= None,track_swaps=False): 
    swaps=np.arange(M_input.shape[0])
=======
def row_echelon(M_input,field, max_col= None): 
>>>>>>> main
    if max_col is None:
        max_col=M_input.shape[1]
    M= Field.flatten_zero(copy.deepcopy(M_input),field)
    #empty matrix
    if M.shape[0]*M.shape[1]==1:
        if M[0][0] == 0:
            return M,[]
        return M,[0]
    
    pivots=[]
    #create one new echelon
    def echelonify(next_pivot_row, col):
<<<<<<< HEAD
        #choose best row to pivot)
=======
        #choose best row to pivot
>>>>>>> main
        if (field.descr in ['Q','R','C']) :
            best_row= next_pivot_row+np.argmax(np.abs(M[next_pivot_row:,col]))
        else:
            non_zero_rows_sub =  next_pivot_row + np.where(   M[next_pivot_row:M.shape[0],col]!= 0  )[0]
            if len(non_zero_rows_sub)==0:
                best_row=next_pivot_row
            else:
                best_row=non_zero_rows_sub[0]
        #swap rows
        if not Field.is_all_zero_elem(M[best_row][col],field):
            rw=np.copy(M[next_pivot_row])
            M[next_pivot_row]=np.copy(M[best_row])
            M[best_row]=rw
            rw=np.copy(M[next_pivot_row,col:])
            pivots.append(col)
<<<<<<< HEAD
            swaps[next_pivot_row],swaps[best_row]=swaps[best_row],swaps[next_pivot_row]
=======
>>>>>>> main
        else: # the column col is null
            return next_pivot_row
        #echelonify the matrix
        non_zero_sub_pivot= next_pivot_row+1+ np.where(   np.transpose(Field.flatten_zero(M[next_pivot_row+1:,col],field)) !=0)[0]
        quotient_values= M[non_zero_sub_pivot,col] /rw[0]
        M[ non_zero_sub_pivot,col:] -= np.matmul( quotient_values.reshape((quotient_values.shape[0],1)) ,rw.reshape((1,rw.shape[0])) )

        
        return next_pivot_row+1

    
    next_pivot_row=0#nb of pivoted rows +1
    for i in range(max_col):#column to pivot
      if next_pivot_row>=M.shape[0]:#all possible rows pivoted
          break
      next_pivot_row=echelonify(next_pivot_row, i)
    #remove some computation errors
    M=Field.flatten_zero(M,field)
    #reduce M
    for i in range(len(pivots)):
        M[i]/= M[i][pivots[i]]
        rw=np.copy(M[i])
        for j in range(i):
            if not(Field.is_all_zero_elem(M[j][pivots[i]],field)):
<<<<<<< HEAD
                M[j] -= rw* M[j][pivots[i]] 
    if track_swaps:
        return np.array(M),pivots,swaps
=======
                M[j] -= M[j][pivots[i]] * rw
>>>>>>> main
    return np.array(M),pivots




#put in column echelon form
<<<<<<< HEAD
def col_echelon(M,field,max_col=None,track_swaps=False):
    row_ech_tr=row_echelon(np.transpose(M),field,max_col=max_col,track_swaps=track_swaps)
    return np.transpose(row_ech_tr[0]),*row_ech_tr[1:]
=======
def col_echelon(M,field,max_col=None):
    row_ech_tr=row_echelon(np.transpose(M),field,max_col=max_col)
    return np.transpose(row_ech_tr[0]),row_ech_tr[1]
>>>>>>> main
     
   
#compute kernel of M
def null_space(M,field):
    if M.shape[1]==0:
        return field.to_Field(np.zeros((0,0),dtype="i"))
    # if field is R or C: SVD
    if field.descr in ['R','C']:
        return scipy.linalg.null_space(M,rcond=epsilon)
    # otherwise column echelon of the augmanted matrix
    M=Field.flatten_zero(M,field)
    aug_mat=Field.flatten_zero(col_echelon(np.concatenate([M, field.to_Field(np.eye(M.shape[1]))]),field)[0],field)
    #column of the kernel base
    zero_col_top=[Field.is_all_zero_mat(aug_mat[:M.shape[0],i],field) for i in range(M.shape[1])]
    return aug_mat[M.shape[0]:,np.array(zero_col_top)]
    
        
#intersection of two families U and V by computing the kernel of 
#( U)
#(-V)
def intersection(U,V,field):
    U=Field.flatten_zero(U,field)
    V=Field.flatten_zero(V,field)
    M=np.concatenate((U,-V),axis=1)
    #empty matrix
    if np.shape(M)[0]*np.shape(M)[1]==0:
        return np.array([]).reshape(np.shape(M))
    u=null_space(M,field)[:np.shape(U)[1]]
    return np.dot(U,u)
    
def extract_basis(M,field):
    if M.shape[0]* M.shape[1]==0:
        return M.reshape((M.shape[0],0))
    if field.descr in['R','C']:
        return scipy.linalg.orth(M,rcond=epsilon)
    col_ech,pivots=col_echelon(M, field)
    return col_ech[:,:len(pivots)].reshape((M.shape[0], len(pivots)))

def complete_basis(M,field):
    if M.shape[1]==0:
        return  field.to_Field(np.eye(M.shape[0]))
    if  field.descr in['Q','R','C']:
        return null_space(np.transpose(M), field)
    #TODO not working when column are swapped
<<<<<<< HEAD
    col_ech,pivots,swaps=col_echelon(np.concatenate((M, field.to_Field(np.eye(M.shape[0]))),axis=1), field,track_swaps=True)
    swaps_inv=np.zeros(len(swaps),dtype="i")
    for i in range(len(swaps)):
        swaps_inv[swaps[i]]=i
    good_cols = [i for i in range(M.shape[0]) if not(Field.is_all_zero_mat(col_ech[:,swaps_inv[M.shape[1]+i]],field))]
    if len(good_cols)==0:
        return field.to_Field(np.zeros((M.shape[0],0),dtype="i"))
    #non_zero_cols=[not(Field.is_all_zero_mat(col_ech[:,M.shape[1]+i],field))  for i in range(M.shape[0])]
    complete_base= field.to_Field(np.eye(M.shape[0]))[:,np.array(good_cols)]
    if extract_basis(np.concatenate([M,complete_base],axis=1), field).shape[1]!=M.shape[0]:
        print(M , "col_ech \n",col_ech,"\n swaps",swaps,"complete_base",complete_base)
        raise ValueError("e3")
    return  complete_base
=======
    col_ech,pivots=col_echelon(np.concatenate((M, field.to_Field(np.eye(M.shape[0]))),axis=1), field)
    non_zero_cols=[not(Field.is_all_zero_mat(col_ech[:,M.shape[1]+i],field))  for i in range(M.shape[0])]
    return  field.to_Field(np.eye(M.shape[0]))[:,np.array(non_zero_cols)]
>>>>>>> main

def sum_subspaces(U,V,field):
    U=Field.flatten_zero(U,field)
    V=Field.flatten_zero(V,field)
    M=np.concatenate((U,V),axis=1)
    return extract_basis(M,field)

# matrix of a projection from dim tot to dim b
def proj(a,b,tot,field,B=None):
    if B==None:
        B= field.to_Field(np.eye(b))
    return np.concatenate((field.to_Field(np.zeros((b,a),dtype="i")),B,field.to_Field(np.zeros((b,tot-a-b),dtype="i")) ),axis=1)

def inj(a,b,tot,field):
    return proj(a,b,tot,field).transpose()

# transform a matrix from row echelon form to diagonal
def ech_to_diag_row(T_input,field):
    T=copy.deepcopy(T_input)
    #P_pivots s.t. T*P_pivot diag
    pivots=[]
    for i in range(min(T.shape)):
        col_piv=i
        while col_piv < T.shape[1] and Field.is_all_zero_elem(T[i][col_piv], field) :
            col_piv+=1
        if col_piv<T.shape[1]:
            pivots.append(col_piv)
        else:
            pivots.append(-1)
    for i in range(min(T.shape)):
        if pivots[i]!=-1:      
            T[i]=T[i]/T[i][pivots[i]]
    for i in range(min(T.shape)):
        if pivots[i]!=-1: 
            for i_2 in range(i):
                T[i_2]= T[i_2] - np.array([T[i_2][pivots[i]]/T[i][pivots[i]]])*T[i]                  
    return T

# transform a matrix from column echelon form to diagonal    
def ech_to_diag_col(T_input,field):
    return np.transpose(ech_to_diag_row(np.transpose(copy.deepcopy(T_input)),field))

# solve Mx=y with M triangular (square)
def solve_triangular(M,y,field):
    if y.shape[1]==0:
        return np.array([]).reshape(M.shape[0],0)
    
    if field.descr in ['R','C']:
        return scipy.linalg.solve_triangular(M,y)
    
    aug_mat=ech_to_diag_row(np.concatenate([M,y],axis=1),field)
    y_ech=aug_mat[:,M.shape[1]:]
    for i in range(M.shape[0]):
        y_ech[i]=y_ech[i]/aug_mat[i][i]
    return y_ech



# compute the inverse image of M restricted to Im(M)\cap K
def inverse_image(M,K,field):
    #column echelon form
    Img_M=col_echelon(M,field)[0]
    #eliminate zero columns
    if field.descr in ['R','C']:
        non_zero_cols=np.where(np.max(np.abs(Img_M),axis=0)>epsilon)[0]
    else:
        non_zero_cols=[i  for i,M_col_i in enumerate(list(Img_M.transpose())) if not(Field.is_all_zero_mat(M_col_i,field)) ]
    Img_M=Img_M[:,non_zero_cols]
    #basis of Im(M)\cap K
    Img_inter=intersection(Img_M,K,field)
    #solve Mx=y for y in Im(M)\cap K
    x=inverse_image_vect(M,Img_inter,field)
    ker = null_space(M,field)#add a basis of ker(M)
    return np.concatenate((x, ker),axis=1)
    
#solve Mx=y
def inverse_image_vect(M,y,field):
    assert(M.shape[0]==y.shape[0])
    #empty matrix
    if M.shape[1]*M.shape[0]==0:
        return np.array([]).reshape(M.shape[1],y.shape[1])
    #column echelon form of the augmented matrix
    Augmented_mat,pivots=row_echelon(np.concatenate((M,y),axis=1),field,max_col=M.shape[1])
    M_ech=Augmented_mat[:,:len(M[0])]
    #no solution
    if len(pivots)>0 and pivots[-1]>= np.shape(M)[1]:
        return np.array([]).reshape(np.shape(M)[0],0)
    #tranform into square invertible triangular matrix and solve
    non_zero_rows = np.array([i for i in range(len(M_ech)) if not(Field.is_all_zero_mat(M_ech[i],field))])
    #M=0
    if len(non_zero_rows)==0:
        assert(Field.is_all_zero_mat(y, field))
        return  field.to_Field(np.eye(M.shape[0]))
    try:
        x_part= solve_triangular(M_ech[non_zero_rows][:,np.array(pivots)],Augmented_mat[non_zero_rows][:,len(M[0]):] ,field)
    except:
        print("erazerze",non_zero_rows,pivots,M,y)
        raise ValueError("e2")
    #reintegrate 0 rows and non pivots
    x=field.to_Field(np.zeros((len(M[0]),len(x_part[0])),dtype="i"))
    for i_p,p in enumerate(pivots):
        x[p]=x_part[i_p]
    return x

#finds a nonzero solution to Mx=y
def inverse_image_vect_from_ech(M_ech,pivots,y,field):
    pivots=np.intersect1d(np.array(pivots), np.array(range(M_ech.shape[1])))
    assert(M_ech.shape[0]==y.shape[0])
    #empty matrix
    if M_ech.shape[1]*M_ech.shape[0] * y.shape[1]==0:
        return np.array([]).reshape(M_ech.shape[1],y.shape[1])
    is_zero_row=np.array([Field.is_all_zero_mat(M_ech[i,:],field) for i in range(M_ech.shape[0])])
    #no solution
    if not(Field.is_all_zero_mat(y[is_zero_row],field)):
        return np.array([]).reshape(np.shape(M_ech)[0],0)
    #extract invertible minor and solve
    M_square,y_square= M_ech[ np.logical_not(is_zero_row)][:,np.array(pivots)], y[ np.logical_not(is_zero_row)]
    x=field.to_Field(np.zeros((M_ech.shape[1], y.shape[1]),dtype="i"))
    for i_p,p in enumerate(pivots):
        x[p]=y_square[i_p]/ M_square[i_p][i_p]
    return x


def spanning_subrep(V,x,Vx):
    Wx=dict(zip(V.vertices, [V.field.to_Field(np.zeros((V.spaces[v], 0),dtype="i")) for v in V.vertices]))
    Wx[x]=Vx
    pile=[x]
    while len(pile)>0:
        for e in V.edges_out[pile.pop()]:
            try:
                Wx[V.edges[e][1]]=Field.flatten_zero(extract_basis(np.dot(V.Ve[e], Wx[V.edges[e][0]]),V.field),V.field)
                assert(Wx[V.edges[e][1]].shape[0]>0 or Wx[V.edges[e][1]].shape[1]==0)
            except:
                sleep(0.1)
                print(np.shape(np.dot(V.Ve[e], Wx[V.edges[e][0]])))
                raise ValueError("erbas")
            pile+=[V.edges[e][1]]
    return Wx
    

def bases_change(V,P):
    assert(V.vertices==list(P.keys()))
    We={}
    for e in V.edges.keys():
        y=np.matmul(V.Ve[e],P[V.edges[e][0]] )
        We[e]= inverse_image_vect(P[V.edges[e][1]], y, V.field)
    return Quiver.Quiver(V.vertices,V.edges,We,V.field,V.grid)
        


def test_subrep_quotient_rep(field_descr):
    field=Field.Field(field_descr)
<<<<<<< HEAD
    V1,_= int_module_in_grid_quiver([5,4], [1,1], [3,2],field)
    V2,_= int_module_in_grid_quiver([5,4], [0,1], [3,3],field)
=======
    V1= int_module_in_grid_quiver([5,4], [1,1], [3,2],field)
    V2= int_module_in_grid_quiver([5,4], [0,1], [3,3],field)
>>>>>>> main
    
    
    V=direct_sum(V1, V2)
    
    
    V.display_graph(label="V")

    
    Wx=spanning_subrep(V, (1,1), np.array([[field.one],[field.one]]))
    
    
    
    W= subrep(V, Wx)
    W.display_graph(label="W")
    
    WT= quotientrep(V, Wx)
    WT.display_graph(label="V/W")

'''
field=Field.Field("Q")
Wx=np.array([1,0,0]).reshape((3,1))*field.one
print(type(complete_basis(Wx, field)[0][0]))
'''


def print_frac_2darray(A):
    for x in range(A.shape[0]):
        for y in range(A.shape[1]):
            print(np.round(float(A[x][y]),2),end="\t")
        print()
        
def compute_quotient_slopes(HN,x_set, vertices):
    for x in x_set:
        for i in range (1,len(HN[x])):
            ind_x= vertices.index(x)
            d_x= HN[x][i][ind_x]-HN[x][i-1][ind_x]
            d_tot= sum( [HN[x][i][y]-HN[x][i-1][y] for y in range (len(HN[x][0]))])
            print(Fraction(int(d_x),int(d_tot)))

def print_grid_HN_type(HN,xmax,x_set):
    for x in x_set:
        print("x=",x)
        for (k,dim_vect) in enumerate(HN[x]):    
            for j in range(xmax[1]-1,-1,-1):
                if np.max(np.int16(dim_vect[j*xmax[0]:(j+1)*xmax[0]]))>0:
                    for i in range (xmax[0]):
                        dim_x=int(dim_vect[i+j*xmax[0]]) 
                        print(dim_x if dim_x!=0 else " " ,end="\t")
                    print()
            print("---")

def ceil_div(x:int, N:int):
    if x%N==0:
        return x//N
    return x//N+1