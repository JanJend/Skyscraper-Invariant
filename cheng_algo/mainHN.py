# -*- coding: utf-8 -*-
import numpy as np
import scipy as scipy
import auxHN as aux
import HNcshrunk as cshrunk
from time import time,sleep
import copy
from cfractions import Fraction
epsilon=1e-10
import warnings
import traceback
import matplotlib.pyplot as plt
import Field
warnings.simplefilter("error")



x_glob=(0,0)

def computeHN_sub(V,x,filtration=False,verbose=False):    # assumes V=<V_x> 
    d=1
    field=V.field

    maps_span= build_spanning_maps(V, x)
    maps_span.pop(x)
    maps_span_l= [ Aj for Aj in maps_span.values() if Aj.shape[0]*Aj.shape[1]!=0] 

        # trivially semistable
    if V.spaces[x] in [0,sum(list(V.spaces.values()))] :
        return[]
    
    v_x=V.spaces[x]
    v_notx=sum(list(V.spaces.values()))-v_x
<<<<<<< HEAD
    #print("v_x=",v_x,"v_notx=",v_notx,"l",list(V.spaces.values()))
    U0=None
    for v_x_small in range(int(min(np.ceil(v_x/float(v_notx)),v_x)),30+2*max(2,v_x+1)):
=======
    U0=None
    for v_x_small in range(int(min(np.ceil(v_x/float(v_notx)),v_x)),max(2,v_x+1)):
>>>>>>> main
        U0=None
        blocks_list= cshrunk.buildblock(maps_span_l,d,v_x_small)
        try:  
            for r in blocks_list:
                blocks,blocks_size=r
                size_X= sum([s[0] for s in blocks_size[:,0]]), sum([s[1] for s in blocks_size[0]])
                M= np.block([[np.eye(min(blocks.shape)),np.zeros((min(blocks.shape), blocks.shape[1]-min(blocks.shape)))],
                             [np.random.randint(1,25*max(size_X) +10, size=(blocks.shape[0]-min(blocks.shape), min(blocks.shape))),  np.zeros( (blocks.shape[0]-min(blocks.shape),  blocks.shape[1]-min(blocks.shape))) ]])
<<<<<<< HEAD
                M=field.to_Field(np.array(M,dtype="i"))
                X= (M,maps_span_l, blocks) 
                A=cshrunk.assemble_block(X)
                if v_x_small>10 and (n & (n - 1)) == 0:
                    print(blocks.shape,"t",aux.extract_basis(A,field).shape[1],aux.extract_basis(A,field).shape[1] *5 / blocks.shape[0]  )          
                try:
                    t=time()
                    U0_temp=cshrunk.wongblockpseudo(X,field=field)
                    t1=time()
                    '''
                    U0p=aux.extract_basis(cshrunk.wongblock(X)[:v_x],field)
                    print(U0_temp.shape,U0p.shape,"e")
=======
                M=field.one*np.array(M,dtype="i")
                X= (M,maps_span_l, blocks)            
                try:
                    t=time()
                    U0_temp=cshrunk.wongblockpseudo(X)
                    t1=time()
                    '''
                    U0p=aux.extract_basis(cshrunk.wongblock(X)[:v_x],field)
>>>>>>> main
                    t2=time()
                    if t1-t>100:
                        print("new", np.round(t1-t,2)," s, old",np.round(t2-t1,2)," s")
                    assert(aux.intersection(U0_temp,U0p,field).shape[1]==max(U0_temp.shape[1],U0p.shape[1]))
                    '''
                except ValueError as e: 
                    #recompute in case we are unlucky with M
                    if str(e) =="nc-rk A<nc-rk B":
                        M= np.block([[np.eye(min(blocks.shape)),np.zeros((min(blocks.shape), blocks.shape[1]-min(blocks.shape)))],
                                    [np.random.randint(1,25*max(size_X) +10, size=(blocks.shape[0]-min(blocks.shape), min(blocks.shape))),  np.zeros( (blocks.shape[0]-min(blocks.shape),  blocks.shape[1]-min(blocks.shape))) ]])
<<<<<<< HEAD
                        M=field.to_Field(np.array(M,dtype="i"))
                        X= (M,maps_span_l, blocks)
                        Ap=cshrunk.assemble_block(X)
                        U0_temp=cshrunk.wongblockpseudo(X,field=field)
=======
                        M=field.one*np.array(M,dtype="i")
                        X= (M,maps_span_l, blocks)
                        U0_temp=cshrunk.wongblockpseudo(X)
>>>>>>> main
                    else: 
                        raise 
                if not( U0 is None or aux.intersection(U0,U0_temp,field).shape[1]==max(U0.shape[1] ,U0_temp.shape[1]) ) :
                    
                    #c_plus
                    span_rep0_temp= aux.spanning_subrep(V, x, U0_temp)
                    subrep0_temp= aux.subrep(V,span_rep0_temp)
                    dims_subrep0_temp=np.array(list(subrep0_temp.spaces.values()))
                    c_plus= (v_notx+v_x)*subrep0_temp.spaces[x]- v_x*np.sum(dims_subrep0_temp)

                    #disc_plus
                    v_notx_s_plus= aux.ceil_div(v_notx*v_x_small,v_x)
                    #v_tot_s_min=  (v_notx*v_x_small)//v_x        +    v_x_small
                    v_x_disc_plus = v_x*(    v_notx_s_plus*subrep0_temp.spaces[x]- v_x_small*(np.sum(dims_subrep0_temp) -subrep0_temp.spaces[x]))
                    disc_plus = v_x_disc_plus/float(v_x_small) 

                     #c_moins
                    span_rep0= aux.spanning_subrep(V, x, U0)
                    subrep0= aux.subrep(V,span_rep0)
                    dims_subrep0=np.array(list(subrep0.spaces.values()))
                    c_moins= (v_notx+v_x)*subrep0.spaces[x]- v_x*np.sum(dims_subrep0)
                    v_notx_s_minus= (v_notx*v_x_small)//v_x
                    v_x_disc_minus = v_x*(    v_notx_s_minus*subrep0.spaces[x]- v_x_small*(np.sum(dims_subrep0) -subrep0.spaces[x]))
<<<<<<< HEAD
=======

                    print("dims Us",U0.shape[1],U0_temp.shape[1], v_x)
>>>>>>> main
                    #compute list of possible ux, u_notx
                    better_U_possible = False
                    for ux in range(U0.shape[1]+1,U0_temp.shape[1]):
                        u_notx_max_plus = (v_notx*ux-c_plus)//v_x 
                        u_notx_max_minus = (v_notx*ux-c_moins)//v_x 
                        u_notx_max = min(u_notx_max_minus,u_notx_max_plus)
                        u_notx_min_plus = aux.ceil_div((v_x-v_x_disc_plus+v_x*ux*(v_notx_s_plus)),v_x_small*v_x)
                        assert( v_x_disc_plus//v_x> 
                                v_notx_s_plus*ux- v_x_small*u_notx_min_plus and 
                                 v_x_disc_plus//v_x <= 
                                    v_notx_s_plus*ux- v_x_small*(u_notx_min_plus-1))
                        u_notx_min_minus = aux.ceil_div((-v_x_disc_minus+v_x*ux*(v_notx_s_minus)),v_x*v_x_small)
                        u_notx_min =max( u_notx_min_minus,  u_notx_min_plus )

                        #print("ux",ux,",",u_notx_min,u_notx_max,",",v_x_small)
                        better_U_possible |=  u_notx_max >= u_notx_min
<<<<<<< HEAD
                    #raise error if bound disc< floor(disc+) not good enough
                    if better_U_possible:
                        raise ValueError("weight approx incorrect")
                    if c_plus<c_moins:
                        U_temp=U0
=======

                    #assert(c_plus>=c_moins)
                    if better_U_possible:
                        pass
                        '''print("U0-=",U0.shape[1],np.sum(dims_subrep0) -U0.shape[1] ,"U0+=",U0_temp.shape[1],
                              np.sum(dims_subrep0_temp) -U0_temp.shape[1],"vx=",v_x,"v_notx=",v_notx,"v_x_s=", v_x_small,"v_notx_s+=",v_notx_s_plus,
                              "c+=",c_plus,"disc+",disc_plus)
                        '''
                    #raise error if bound disc< floor(disc+) not good enough
                    if better_U_possible:
                    #if  c_plus_app *v_x_small < v_x_disc_plus  :    
                        if v_x>2:
                            pass
                            #print(size_X,v_x,v_notx)  
                            #print("c-=",c_moins,end=" | ")
                            #print("c+=",c_plus," ","disc+=", disc_plus, end =" ")
                            #print("c+_app=", c_plus_app,end=" " )
                            #print("disc-=", (v_x/float(v_x_small))  *(     v_tot_s_min*subrep0.spaces[x]- v_x_small*np.sum(dims_subrep0)))
                            #print("U0-")
                            #aux.print_frac_2darray(U0) 
                            #print("U0+")
                            #aux.print_frac_2darray(U0_temp)
                        raise ValueError("weight approx incorrect")
                    if c_plus<c_moins:
                        U_temp=U0
                    if v_x>2:
                        pass#print("c+=",c_plus," /\ ","disc+=", disc_plus,"c+_app=",c_plus_app, "v_x_small=",v_x_small, end =" ")
                    #print("ux=",U0_temp.shape[1],"v_x_small",v_x_small, "v_s_notx=", (v_notx*v_x_small)//v_x    +1     , "v_x=",v_x,"v_notx=",v_notx)
>>>>>>> main
                U0=U0_temp
            break
        except ValueError as e:
            if str(e) =="weight approx incorrect":
<<<<<<< HEAD
                if v_x>2 and v_x_small>0 and verbose:
                    print("weigt err:",v_x_small,v_x)
            elif str(e)== "nc-rk A<nc-rk B":
                print("blow-up err:",v_x_small,v_x)
                if v_x_small>v_x:
                    pass
=======
                if v_x>2 and v_x_small>0:
                    print("weigt err:",v_x_small,v_x)
            elif str(e)== "nc-rk A<nc-rk B":
                print("blow-up err:",v_x_small,v_x)
>>>>>>> main
            else: 
                raise 

    #print(aux.intersection(Ublock,U,V.field).shape[1]==U.shape[1])


    U0=Field.flatten_zero(U0, field)
    if v_x>3 and v_x_small>1:#verbose:
        pass
        #print("U0=")
        #print(aux.print_frac_2darray(U0))
    if U0.shape[1]==V.spaces[x]:
        return[]
    span_rep= aux.spanning_subrep(V, x, U0)
    subrep= aux.subrep(V,span_rep)
    dims_subrep=np.array(list(subrep.spaces.values()))
    if v_x>6 and v_x_small>3:
        print("c=",(v_notx+v_x)*subrep.spaces[x]- v_x*np.sum(dims_subrep))
    l=[]
    if verbose:
        print("srep",np.sum(dims_subrep))
    if not(subrep.is_zero()):
        l=computeHN_sub(subrep,x,filtration,verbose)+[dims_subrep]
        quotrep= aux.quotientrep(V, span_rep)
        if verbose:
            print("qrep",sum(quotrep.spaces.values()))
        if not(quotrep.is_zero()):
<<<<<<< HEAD
            try:
                l=l+ [dims_subrep+ dims_quot for dims_quot in computeHN_sub(quotrep,x,filtration,verbose)]
            except ValueError as e:
                V.display_graph("V")
                subrep.display_graph("U")
                quotrep.display_graph("V/U")
                print('V', V.spaces.values(),'subrep', subrep.spaces.values(), 'quotrep', quotrep.spaces.values())
                raise
=======
            l=l+ [dims_subrep+ dims_quot for dims_quot in computeHN_sub(quotrep,x,filtration,verbose)]
>>>>>>> main
    return l
    
def computeHN(V,x_set=False,filtration=False,verbose=False): 
    if x_set==False:
        x_set=V.vertices
    skyscraper={}
    for x in x_set:
        if V.spaces[x]==0:
            skyscraper[x]= [np.zeros(len(V.vertices),dtype="i"),np.array(list(V.spaces.values()))]
        else:
            span_maps= aux.spanning_subrep(V, x,  V.field.to_Field(np.eye(V.spaces[x])))
            span_subrep= aux.subrep(V,span_maps)
<<<<<<< HEAD
            try:
                l=computeHN_sub(span_subrep,x,filtration,verbose)
            except ValueError as e:
                V.display_graph("V2")
                span_subrep.display_graph("U2")
                raise
=======
            l=computeHN_sub(span_subrep,x,filtration,verbose)
>>>>>>> main
            V.display_graph("V",verbose)
            skyscraper[x]= [np.zeros(len(V.vertices),dtype="i")]+l+[np.array(list(span_subrep.spaces.values()))]
            if sum(list(span_subrep.spaces.values()))<sum(list(V.spaces.values())):
                span_subrep.display_graph("<V_{"+str(x)+"}>",verbose)
                skyscraper[x]+= [np.array(list(V.spaces.values()))]
    return skyscraper



def build_spanning_maps(V,x):
    field=V.field
    maps_span={x: V.field.to_Field(np.eye(V.spaces[x]))}
    file=[x]
    while len(file)>0:
        for e in V.edges_out[file.pop()]:
            if  V.edges[e][1]  in maps_span.keys():
                #check commutativity
                try:
                    assert(Field.is_all_zero_mat(maps_span[V.edges[e][1]]-np.dot( V.Ve[e], maps_span[V.edges[e][0]]), field))
                except:
<<<<<<< HEAD
                    V.display_graph("V")
                    print(x)
                    print("r",V.edges[e],maps_span[V.edges[e][1]],np.dot( V.Ve[e], maps_span[V.edges[e][0]]), field)
                    raise
            else: 
                file.append(V.edges[e][1])
                try:
                    maps_span[V.edges[e][1]]=np.dot( V.Ve[e], maps_span[V.edges[e][0]])
                except:
                    print("e",V.Ve[e], maps_span[V.edges[e][0]])
                    raise ValueError("error in map "+str(e)+" from "+str(V.edges[e][0])+" to "+str(V.edges[e][1]))
=======
                    print("r",V.edges[e],maps_span[V.edges[e][1]],np.dot( V.Ve[e], maps_span[V.edges[e][0]]), field)
            else: 
                file.append(V.edges[e][1])
                maps_span[V.edges[e][1]]=np.dot( V.Ve[e], maps_span[V.edges[e][0]])
>>>>>>> main
    return maps_span


def random_change_bases(V):
    P={}
    for x in V.vertices:
        if V.spaces[x]==0:
            P[x]=V.field.to_Field(np.ones((0, 0)))
        else:
            p= np.random.randint(1,20,(V.spaces[x],V.spaces[x]))
<<<<<<< HEAD
            while (V.field.descr in ['Q','R','C'] and
                    not scipy.linalg.orth(p,rcond=epsilon).shape[1]==V.spaces[x]) or \
                    (V.field.descr not in ['Q','R','C'] and aux.extract_basis(V.field.to_Field(p),V.field).shape[1]!=V.spaces[x]):
                p= np.random.randint(1,20,(V.spaces[x],V.spaces[x]))
            P[x]=V.field.to_Field(p)
    return aux.bases_change(V,P)
    
def test_skyscraper(grid_size,n_int,verbose=False,field=Field.Field("Q") ):
    
    
=======
            porth=scipy.linalg.orth(p,rcond=epsilon)
            while porth.shape[1]!=porth.shape[0]:
                p= np.random.randint(1,20,(V.spaces[x],V.spaces[x]))
                porth=scipy.linalg.orth(p,rcond=epsilon)
            P[x]=p*V.field.one
    return aux.bases_change(V,P)
    
def test_skyscraper(grid_size,n_int,verbose=False ):
    
    field=Field.Field("Q")
>>>>>>> main
    V,_=aux.int_module_in_grid_quiver(grid_size, [(1,1)], [(0,0)],field)
    list_V=[]
    if n_int>2:
        V=aux.grid_indec(grid_size,field)
        n_int-=2
        list_V=[copy.deepcopy(V)]
    
    for k in range(n_int):
        s=np.random.randint(np.zeros(len(grid_size)),1#(2*np.array(grid_size))//4
                            ,2)
        t0=np.random.randint(s,np.array(grid_size),2)
        t1=np.random.randint(s,np.array(grid_size),2)
        if verbose:
            print("source,taget: ",s,t0,t1)
        V1,supp1= aux.int_module_in_grid_quiver(grid_size, [tuple(s)], [tuple(t0),tuple(t1)],field)
        #V1,supp1= aux.int_module_in_grid_quiver(grid_size, [(min(s0),min(s1))], [(max(s0),max(s1))],field)
        V=aux.direct_sum(V, V1)
        list_V.append(V1)
    #V.display_graph("V")
    
    W=random_change_bases(V)

    #W.display_graph("W")
    x_set=False
    #x_set=[(0,0)]
    #np.random.seed(2)
    HN1= compute_HN_from_support(list_V,x_set=x_set)
    HN2=computeHN(W,x_set=x_set   ,verbose=verbose)
    #HN2bis= computeHN(W,x_set=x_set   ,verbose=False)
    


    
    assert(list(HN1.keys())==list(HN2.keys()))
    success=True
    for x0 in HN1.keys():
        try:
            assert(len(HN1[x0])==len(HN2[x0]) and all([Field.is_all_zero_mat(HN1[x0][i]-HN2[x0][i],V.field) for i in range(len(HN1[x0]))] ))
        except:
            print(x0)
            print("computing from supports:")
            print(HN1[x0])
            print("computing from HN:")
            print(HN2[x0])
            #print(len(HN1[x0])==len(HN2[x0]),[aux.is_all_zero_mat(HN1[x0][i]-HN2[x0][i],V.field) for i in range(len(HN1[x0]))])
            raise ValueError("i")
            success=False
    return success

    

    
def x_span_set(vertices,x):
    return dict(zip(vertices,[all(np.array(v)>=np.array(x)) for v in vertices] ))

def intersect_set(S1,S2):
    assert(list(S1.keys())==list(S2.keys()))
    return dict(zip(S1.keys(),[S1[v] and S2[v] for v in list(S1.keys())] ))

def sum_HN(HN1,HN2,vertices):
    null_dim_vect= np.zeros(len(vertices),dtype="i")
    cards=set(list(HN1.keys())+ list(HN2.keys()))
    HN= dict(zip(cards ,  [np.copy(null_dim_vect) for _ in cards])  )
    for card in list(HN1.keys()):
        HN[card]+=  HN1[card]
    for card in list(HN2.keys()):
        HN[card]+=  HN2[card]
    return HN




def compute_HN_from_support(V,x_set=False):
    vertices=V[0].vertices
    if not(x_set):
        x_set=vertices
    skyscraper={}
    for x in x_set:
        totdims=np.sum([ np.array( list(W.spaces.values()),dtype=np.int32) for W in V ],axis=0)
        HN={}
        #interval modules
        for W in V:
            if W.spaces[x]>0:
                HN= sum_HN( HN, HN_everywhere_ss(W,x) ,vertices)      
        quotient_dims = [a[1] for a in  sorted(HN.items())]
        filtr_dims=copy.copy(quotient_dims)
        for i in range(len(filtr_dims)-1):
            filtr_dims[i+1]+=filtr_dims[i] 
        skyscraper[x]=[np.zeros_like(totdims)]+filtr_dims
        if np.any( totdims-skyscraper[x][-1]>0):
            skyscraper[x]+=[totdims]
    return skyscraper

def HN_everywhere_ss(V,x):
        support=dict(zip(V.vertices,[V.spaces[v]!=0 for v in V.vertices]))
        assert( x in support)    
        S_cap_x=intersect_set(support, x_span_set(V.vertices, x))
        dim_intersect = np.array([  int(S_cap_x[y])*V.spaces[y ] for  y in V.vertices],dtype="i")
        card= np.sum(dim_intersect)/float(V.spaces[x])
        return {card:  dim_intersect} 

<<<<<<< HEAD

field=Field.Field("F_2")
n=5
xmax=[n,n]
t= time()
success_fail=[0,0]
#start_test=np.random.randint(0,10000)
start_test = 154540+60#
n_test=10
print("random seed=",start_test)
for k in range(0,n_test):
    np.random.seed(start_test+k)
    print("k=",k)

    success_fail[test_skyscraper(xmax,n,False,field)]+=1

print(np.round((time()-t)/n_test,2),"s")
print("success: ",int(100* success_fail[1]/sum(success_fail)),"%")


'''

for k in range(10):
    print("k=",k)
    x_set=[(0,0)]
    V= aux.random_grid_indec(xmax, 4,3, 2)
=======
'''
t= time()
success_fail=[0,0]
start_test=np.random.randint(0,10000)#154540+60#
n_test=10
print("random seed=",start_test)
for k in range(n_test):
    np.random.seed(start_test+k)
    print("k=",k)
    n=6
    success_fail[test_skyscraper([n,n],n,False)]+=1

print(np.round((time()-t)/n_test,2),"s")
print("success: ",int(100* success_fail[1]/sum(success_fail)),"%")
'''

field=Field.Field("Q")
xmax=[10,10]

for k in range(10):
    print("k=",k)
    x_set=[(0,0)]
    V= aux.random_grid_indec(xmax, 4,4, 2)
>>>>>>> main
    #x_set=[0]
    #V= aux.star_quiver(5, 8, 3)
    
    V=random_change_bases(V)
    #V.display_graph("V")
    
    
    HN=computeHN(V,x_set)
    #aux.print_grid_HN_type(HN, xmax,  x_set)
    aux.compute_quotient_slopes(HN, x_set, V.vertices)
    #print(len(HN[(0,0)]))
    print("len=",len(HN[x_set[0]]))
<<<<<<< HEAD
=======
'''
>>>>>>> main

V1,supp1= aux.int_module_in_grid_quiver( xmax, [(0,0)],[(0,2),(1,1)],field)
V2,supp2= aux.int_module_in_grid_quiver( xmax, [(0,0)],[(2,1),(1,2)],field)
#V3,supp3= aux.int_module_in_grid_quiver( xmax, [(1,1)],[(0,0)],field)
V3=aux.grid_indec(xmax,field)

V=aux.direct_sum(aux.direct_sum(V1, V2), V3)
#V.display_graph("V before change of basis")
V=random_change_bases(V)



HN1= compute_HN_from_support( [ V1,V2,V3],False)
HN2=computeHN(V)
print("computing from HN:")
print(HN2[(0,0)])
print("computing from supports:")
print(HN1[(0,0)])
<<<<<<< HEAD
'''
=======
'''


>>>>>>> main
