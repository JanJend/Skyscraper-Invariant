import numpy as np
import scipy as sc
import copy as copy
import auxHN as aux
from time import time,sleep
import math
import Field
import Quiver
epsilon=1e-10





<<<<<<< HEAD
def wongblock(X,field=Field.Field("Q")):
    
    M,V,blocks=X 
    X_ass= assemble_block(X)   
    vx= V[blocks[0][0] ].shape[1]
    #i=0
    ind=0
    Wold,W = (None,field.to_Field(np.zeros((X_ass.shape[0],0),dtype="i")))
    U= ker_block_Q(X,field)
    while ind==0 or W.shape[1]>Wold.shape[1]:
        ind+=1
        Wold=W
        W_list=[]
        for i in range(len(V)):
            W_list.append( aux.extract_basis(np.block([ np.matmul(V[i], U[vx*y:vx*(y+1)]) for y in range(len(U)//vx)  ]),field))
        W_list= [ W_list [i%len(V)] for i in range(blocks.shape[0])]
        W= Field.build_block_diag_l(W_list ,field)
        WcapImA= aux.intersection(W,aux.extract_basis(X_ass,field),field)
        if WcapImA.shape[1]!= W.shape[1]:
            raise ValueError("nc-rk A<nc-rk B")
        U=aux.inverse_image(X_ass, W, field)    
        
    
    '''
    kerA=ker_block_Q(X,field)
    vx= V[blocks[0][0] ].shape[1]
    
=======
def wongblock(X):
    field=Field.Field("Q")
    
    M,V,blocks=X    
    
    kerA= ker_block_Q(X)
    vx= V[blocks[0][0] ].shape[1]
>>>>>>> main
    W1_list=[]
    for i in range(len(V)):
        W1_list.append( aux.extract_basis(np.block([ np.matmul(V[i], kerA[vx*y:vx*(y+1)]) for y in range(len(kerA)//vx)  ]),field))
    W1_list= [ W1_list [i%len(V)] for i in range(blocks.shape[0])]
    W1= Field.build_block_diag_l(W1_list ,field)
    X_ass= assemble_block(X)
    W1capImA= aux.intersection(W1,aux.extract_basis(X_ass,field),field)
    if W1capImA.shape[1]!= W1.shape[1]:
        raise ValueError("nc-rk A<nc-rk B")
    U=aux.inverse_image(X_ass, W1, field)    

    W2_list=[]
    for i in range(len(V)):
        W2_list.append( aux.extract_basis(np.block([ np.matmul(V[i], U[vx*y:vx*(y+1)]) for y in range(len(U)//vx)  ]),field))
    W2_list= [ W2_list [i%len(V)] for i in range(blocks.shape[0])]
    W2= Field.build_block_diag_l(W2_list ,field)
    assert (W2.shape[1]== W1.shape[1])
<<<<<<< HEAD
    '''
    
    U0=np.concatenate([U[vx*i:vx*(i+1)] for i in range(0,len(U)//vx) ], axis=1)
    U0 =aux.extract_basis(U0,field)
    return U0
=======
    return U
>>>>>>> main











<<<<<<< HEAD
def wongblockpseudo(X,field=Field.Field("Q"),is_M_ech=True):

=======
def wongblockpseudo(X,is_M_ech=True):
    field=Field.Field("Q")
>>>>>>> main
    M,V,blocks=X   
    A= assemble_block(X)
    vx= V[blocks[0][0] ].shape[1]

    

    t0=time()
    #column echelon M
    M_ech , PM= M , field.to_Field(np.eye(M.shape[1]))
    if not(is_M_ech):
        aug_mat_M=aux.ech_to_diag_col(aux.col_echelon(np.concatenate([M, field.to_Field(np.eye(M.shape[1]))]),field)[0],field)
        M_ech , PM = aug_mat_M[:M.shape[0]],aug_mat_M[M.shape[0]:]
    
            #row echelon V
    V_ech_row,PV_row=[],[]
    for Vy in V:   
        aug_mat=aux.row_echelon(np.concatenate([Vy,  field.to_Field(np.eye(Vy.shape[0]))],axis=1),field)[0]
        Vy_ech,Py = aug_mat[:,:Vy.shape[1]],aug_mat[:,Vy.shape[1]:]
        V_ech_row.append(Vy_ech)
        PV_row.append(Py)
<<<<<<< HEAD
            
    kerA = ker_block_Q(X,field)
=======

    
    '''
    # column echelon V
    V_ech_col,PV_col=[],[]
    zero_col_top_V=[]
    for Vy in V:   
        aug_mat=aux.col_echelon(np.concatenate([Vy, field.to_Field(np.eye(Vy.shape[1]))]),field)[0]
        Vy_ech,Py = aug_mat[:Vy.shape[0]],aug_mat[Vy.shape[0]:]
        V_ech_col.append(Vy_ech)
        PV_col.append(Py)
        zero_col_top_V.append(np.where([Field.is_all_zero_mat(Vy_ech[:,i],field) for i in range(Vy_ech.shape[1])])[0])

    #aug_mat_col,pivots_col= aux.col_echelon(np.concatenate([A, field.to_Field(np.eye(A.shape[1]))]), field)
    aug_mat_col = col_ech_block(X,V_ech_col,PV_col,PM)
    A_col_ech,P_col_ech =aug_mat_col[:A.shape[0]],aug_mat_col[A.shape[0]:]
    #is_zero_col=[Field.is_all_zero_mat(A_col_ech[:,i],field) for i in range(A_col_ech.shape[1])]
    assert(Field.is_all_zero_mat( A_col_ech-np.matmul(A,P_col_ech),field))
    '''

    
    kerA = ker_block_Q(X)
>>>>>>> main

    t1=time()



<<<<<<< HEAD
    aug_mat_row,pivots_row= row_echelon_block(X,V_ech_row,PV_row,field) 
=======
    aug_mat_row,pivots_row= row_echelon_block(X,V_ech_row,PV_row) 
>>>>>>> main
    A_row_ech,P_row_ech= aug_mat_row[:,:A.shape[1]],aug_mat_row[:,A.shape[1]:]

    is_zero_row = np.all(A_row_ech==0,axis=1)


    t2=time()

    
    

    #kerA2 = P_col_ech2[:,np.array(is_zero_col2)]

    #assert( aux.intersection(kerA, kerA2,field).shape[1]== max(kerA.shape[1],kerA2.shape[1] ) )
    #imA = A_col_ech[:, np.logical_not(is_zero_col)  ] 
    
<<<<<<< HEAD
    def W_from_U(U):
        W_list=[]
        for i in range(len(V)):
            W_list.append( aux.extract_basis(np.block([ np.matmul(V[i], U[vx*y:vx*(y+1)]) for y in range(len(U)//vx)  ]),field))
        W_list=[W_list [i%len(V)].copy() for i in range(blocks.shape[0])]
        return sum([wi.shape[1] for wi in W_list]),Field.build_block_diag_l(W_list ,field)
    
    #W1
    U=kerA
    W_shape,W= W_from_U(kerA)
    '''
    U0= aux.extract_basis(kerA[:vx],field)
=======
    #W0
    U0= aux.extract_basis(kerA[:vx],field)

    
>>>>>>> main
    W_list=[]
    for i in range(len(V)):
        W_list.append( aux.extract_basis(np.block([ np.matmul(V[i], kerA[vx*y:vx*(y+1)]) for y in range(len(kerA)//vx)  ]),field))
    W_list= [ W_list [i%len(V)] for i in range(blocks.shape[0])]
    W_shape1=sum([ W_list [i%len(V)].shape[1] for i in range(blocks.shape[0])])
<<<<<<< HEAD
    '''
    
    is_W_stab= (W_shape==0)
    ind= 1
    while not(is_W_stab) and ind <= min(A.shape):
        #W= Field.build_block_diag_l(W_list ,field)#W_i
=======
    
    
    is_W_stab= (W_shape1==0)
    ind= 1

    while not(is_W_stab) and i <= min(A.shape):
        W= Field.build_block_diag_l(W_list ,field)#W_i
>>>>>>> main

        PW= np.matmul(P_row_ech, W)
        if not(Field.is_all_zero_mat(PW[is_zero_row],field)):#check W\subset Im A
            raise ValueError("nc-rk A<nc-rk B")

        U= np.concatenate([aux.inverse_image_vect_from_ech(A_row_ech,pivots_row, PW ,field), kerA],axis=1)    
<<<<<<< HEAD
        
        W_shape_old=W_shape
        W_shape,W= W_from_U(U)
        is_W_stab |= (W_shape== W_shape_old)
        ind+=1
        '''
        U0=np.concatenate([U[vx*i:vx*(i+1)] for i in range(1,len(U)//vx) ], axis=1)
        U0 =aux.extract_basis(U0,field)
        
      #Wi+1        
=======
        U0=np.concatenate([U[vx*i:vx*(i+1)] for i in range(1,len(U)//vx) ], axis=1)
        U0 =aux.extract_basis(U0,field)
        #check U= U0 \otimes C^n
        '''
        for i in range(1,len(U)//vx):
            U0check=aux.extract_basis(U[vx*i:vx*(i+1)],field)
            assert(aux.intersection(U0,U0check,field).shape[1]==max(U0.shape[1],U0check.shape[1]))
        '''
        #Wi+1        
>>>>>>> main
        W2_list= [aux.extract_basis( np.matmul(Vy, U0),field) for Vy in V]
        W2_shape1=sum([ W2_list [i%len(V)].shape[1] for i in range(blocks.shape[0])])
        is_W_stab |= (W2_shape1== W.shape[1])
        W_list=W2_list
        ind+=1
<<<<<<< HEAD
        '''
=======
    if False:
        W= aux.build_block_diag_l(W_list ,field)#W_i
    
        PW= np.matmul(P_row_ech, W)
        if not(Field.is_all_zero_mat(PW[is_zero_row],field)):#check W\subset Im A
            raise ValueError("nc-rk A<nc-rk B")
    
        U= np.concatenate([aux.inverse_image_vect_from_ech(A_row_ech,pivots_row, PW ,field), kerA],axis=1)    
        U0=aux.extract_basis(U[0:vx],field)
        #check U= U0 \otimes C^n
        for i in range(1,len(U)//vx):
            U0check=aux.extract_basis(U[vx*i:vx*(i+1)],field)
            assert(aux.intersection(U0,U0check,field).shape[1]==max(U0.shape[1],U0check.shape[1]))
    
    
        W2_list=[]
        for i in range(len(V)):
            W2_list.append( aux.extract_basis(np.block([ np.matmul(V[i], U[vx*y:vx*(y+1)]) for y in range(len(U)//vx)  ]),field))
        W2_list= [ W2_list [i%len(V)] for i in range(blocks.shape[0])]
        W2= aux.build_block_diag_l(W2_list ,field)
        assert (W2.shape[1]== W.shape[1])
>>>>>>> main

    t3=time()
    if t3-t0>100:
        print("ech col", np.round(t1-t0,2), "ech row", np.round(t2-t1,2), "comp", np.round(t3-t2,2))
        if ind>2:
            print("i=",ind)
<<<<<<< HEAD
    if U.shape[1]==0:
        return U[:vx]
    #check U= U0 \otimes C^n
    U0=np.concatenate([U[vx*i:vx*(i+1)] for i in range(0,len(U)//vx) ], axis=1)
    for i in range(0,len(U)//vx):
        U0check=aux.extract_basis(U[vx*i:vx*(i+1)],field)
        assert(aux.intersection(U0,U0check,field).shape[1]==max(U0.shape[1],U0check.shape[1]))
    U0 =aux.extract_basis(U0,field)   
=======
        
>>>>>>> main
    return U0










def assemble_block(X):
    M,V,blocks=X
    if blocks.shape[0]*blocks.shape[1]==0:
        return np.zeros_like(blocks)
    blocks_map= []
    for i in range (blocks.shape[0]):
        blocks_map_i=[]
        for j in range(blocks.shape[1]):
<<<<<<< HEAD
            blocks_map_i.append(  V[blocks[i][j]]* M[i][j])
=======
            blocks_map_i.append( M[i][j]* V[blocks[i][j]])
>>>>>>> main
        blocks_map.append(blocks_map_i)

    return np.block(blocks_map)




def E(n,m,i,j,field):
    E=field.to_Field(np.zeros((n,m),dtype="i"))
    E[i][j]=field.one
    return E

def buildA(V,d,field):
    vx=np.shape(V[0])[1]
    vtot = sum( [np.shape(Vy)[0] for Vy in V])
    Al=[]
    zero_like_V= [field.to_Field(np.zeros(Vy.shape,dtype="i")) for Vy in V]
    for y in range(len(V)):
        Aylist= copy.copy(zero_like_V)
        Aylist[y]=V[y]
        Al.append(np.concatenate(Aylist))
    Al2=[]
    for A in Al:
        for i in range(vx ):
            for j in range (vtot ):
                Al2.append( np.kron( E(vx,vtot,i,j,field),A)) 
    Ad=[]
    for A in Al2:
        for j in range(d):
            for i in range (d):
                Ad.append( np.kron(E(d,d,i,j,field),A)) 
    return Ad

def buildblock(V,d,v_x_small=1):
    
<<<<<<< HEAD

=======
    field=Field.Field("Q")
>>>>>>> main
    v_x=np.shape(V[0])[1]
    v_notx = sum( [np.shape(Vy)[0] for Vy in V])


    
    n_blocks_j= [ (v_notx*d*v_x_small)//v_x] 
    if  (v_notx*d*v_x_small)%v_x!=0:
        n_blocks_j.append((v_notx*d*v_x_small)//v_x+1)

    res=[]
    for nbj in n_blocks_j:
        blocks=np.zeros((v_x_small*len(V)*d,nbj),dtype="i")
        blocks_sizes=np.zeros((v_x_small*len(V)*d,nbj),dtype="i,i")
        for i in range(v_x_small*len(V)*d):
            y=i% len(V)
            for j in range(nbj):
                blocks[i][j],blocks_sizes[i][j]=y,V[y].shape
        res.append((blocks,blocks_sizes))
    return res




def pseudo_inverse(A):
    t=time()
    field=Field.Field("Q")
    kerA = aux.null_space(A,field)
    U=aux.complete_basis(kerA,field)
    P=np.concatenate((U,kerA),axis=1)
    ImA=np.matmul(A,U)
    Up=aux.complete_basis(ImA,field)
    Q=np.concatenate((ImA,Up),axis=1)
    Ap=np.transpose(aux.inverse_image_vect(np.transpose(Q),np.transpose(P),field))
    return Ap


<<<<<<< HEAD
def ker_block_Q(X,field=Field.Field("Q")):
    (M,V,blocks)=X
    #print("blocks =\n",blocks.shape)
    #print("V=\n",V)
    #print("M=\n",M)
    #print("A=","\n",assemble_block(X))
=======
def ker_block_Q(X):
    (M,V,blocks)=X
    #print("A=","\n",assemble_block(X))
    #print("M=\n",M)
    field=Field.Field("Q")
>>>>>>> main
    shape_X= ( sum([ V[blocks[i][0] ].shape[0] for i in range(blocks.shape[0]) ])
        ,sum([V[blocks[0][j]].shape[1] for j in range(blocks.shape[1])])  )
        #column echelon M
    aug_mat=aux.ech_to_diag_col(aux.col_echelon(np.concatenate([M, field.to_Field(np.eye(M.shape[1]))]),field)[0],field)
    M_ech , PM = aug_mat[:M.shape[0]],aug_mat[M.shape[0]:]
    zero_col_top_M=np.where([Field.is_all_zero_mat(M_ech[:,i],field) for i in range(M.shape[1])])[0]
    #print("M_ech=\n",M_ech)
        #column echelon Vy
        # column echelon V
    V_ech,PV=[],[]
    zero_col_top_V=[]
    for Vy in V:   
        aug_mat=aux.col_echelon(np.concatenate([Vy, field.to_Field(np.eye(Vy.shape[1]))]),field)[0]
        Vy_ech,Py = aug_mat[:Vy.shape[0]],aug_mat[Vy.shape[0]:]
        V_ech.append(Vy_ech)
        PV.append(Py)
        zero_col_top_V.append(np.where([Field.is_all_zero_mat(Vy_ech[:,i],field) for i in range(Vy_ech.shape[1])])[0])




    #ViPj=[np.matmul(V_ech[i% ], PV[j]) for j in range()  ]
        # column echelon X
    #X_ech=(M_ech,V_ech,blocks)
    blocksP= np.array( [np.remainder(np.arange(blocks.shape[1]),len(V)) for _ in range(blocks.shape[1])])
    P = (PM, PV, blocksP)
    #print("X_ech=\n", np.matmul(assemble_block(X), assemble_block(P)))
        #kernel from block M
    ker = field.to_Field(np.zeros((shape_X[1],0),dtype="i"))
    if len(zero_col_top_M)>0:
        ker= assemble_block((PM[:,zero_col_top_M], PV, blocksP[:,zero_col_top_M] ))
        # matrix of interdependence
    cols_B_mod=[ y for y in range (len(V)) if len(zero_col_top_V[y])>0 ] 
    cols_B=[ y for y in range (blocks.shape[1]) if len(zero_col_top_V[y%len(V)])>0 ] #col ker X
    if len(cols_B)==0:
        return ker
    QV=[PV[y][:,zero_col_top_V[y]] for y in cols_B_mod ]
    QM=PM[:,cols_B]

    blocksQ=np.empty(QM.shape,dtype=object)
    for i in range(QM.shape[0]):
        for j in range(QM.shape[1]):
            blocksQ[i][j]= j%len(cols_B_mod)


    VxPy=[]
    for y in range(len(cols_B_mod)):
        cols=zero_col_top_V[cols_B_mod[y]]
        for x in range(len(V)):   
            VxPy.append( np.matmul(V[x], QV[y] ))
    
    B=M_ech[M_ech.shape[1]:,cols_B]

    
    blocksB=np.empty((blocks.shape[0]-blocks.shape[1],len(cols_B)),dtype=object)
    for j in range(len(cols_B)):
        for i in range(blocksB.shape[0]):
            (x,y)=((i+blocks.shape[1])%len(V),j%len(cols_B_mod))
            blocksB[i][j]= len(V)*y+x
    
    B_assembled=assemble_block((B,VxPy,blocksB))
    shape_1_B_ass=sum([len(zero_col_top_V[y%len(V)]) for y in cols_B ])
    if blocksB.shape[0]==0:
        B_assembled = field.to_Field(np.zeros((0,shape_1_B_ass),dtype="i"))
    assert(B_assembled.shape[1]== shape_1_B_ass)
    B_augm=aux.col_echelon(np.concatenate((B_assembled, assemble_block((QM,QV,blocksQ)))),field)[0]
    B_ech , Q2 = B_augm[:B_assembled.shape[0]],B_augm[B_assembled.shape[0]:]
    zero_col_top_B=np.where([Field.is_all_zero_mat(B_ech[:,i],field) for i in range(B_assembled.shape[1])])[0]
    #print(Q2[:,zero_col_top_B])
    return np.concatenate([ker,Q2[:,zero_col_top_B]],axis=1)


#column echelon X from column-echeloned M  and V
def col_ech_block(X_input,V_ech,PV,PM):
    X=copy.deepcopy(X_input)
    (M,V,blocks)=X
    field=Field.Field("Q")
    shape_X= ( sum([ V[blocks[i][0] ].shape[0] for i in range(blocks.shape[0]) ])
        ,sum([V[blocks[0][j]].shape[1] for j in range(blocks.shape[1])])  )
        # column echelon the diagonal part of X
    VxPy_full=[]
    for x in range(len(V)):
        for y in range(len(V)):
            VxPy_full.append( np.matmul (V[x], PV[y])) 
    for i in range(blocks.shape[0]):
        for j in range(blocks.shape[1]):
            blocks[i][j]= (i%len(V))*len(V)+ j%len(V)
    X=(M, VxPy_full,blocks) 
        # augmented matrix of X
    blocksP= np.transpose(np.array( [np.remainder(np.arange(blocks.shape[1]),len(V)) for _ in range(blocks.shape[1])]))
    P = (PM, PV, blocksP)
        #null columns of Vy
    zero_col_top_V=[]
    for y in range(len(V)):   
        zero_col_top_V.append(np.where([Field.is_all_zero_mat(V_ech[y][:,i],field) for i in range(V_ech[y].shape[1])])[0])
    cols_B_mod=[ y for y in range (len(V)) if len(zero_col_top_V[y])>0 ] 
    cols_B=[ y for y in range (blocks.shape[1]) if len(zero_col_top_V[y%len(V)])>0 ] #col ker X
        #if A already reduced
    if len(cols_B_mod)==0 or blocks.shape[0]<=blocks.shape[1]:
        return np.concatenate([assemble_block(X), assemble_block(P)])
    # matrix of interdependence
    VxPy_red=[ W[:,zero_col_top_V[k%len(V) ]] for k,W in enumerate(VxPy_full)    ]
    B=M[M.shape[1]:,cols_B]  
    blocksB=np.empty((blocks.shape[0]-blocks.shape[1],len(cols_B)),dtype=object)
    for j in range(len(cols_B)):
        for i in range(blocksB.shape[0]):
            (x,y)=((i+blocks.shape[1])%len(V),cols_B[j]%len(V))
            blocksB[i][j]= len(V)*x+y
    B_ass= assemble_block((B,VxPy_red,blocksB))
    #augmented matrix of B
    QV=[PV[y][:,zero_col_top_V[y]] for y in cols_B_mod ]
    QM=PM[:,cols_B]
    blocksQ=np.empty(QM.shape,dtype=object)
    for i in range(QM.shape[0]):
        for j in range(QM.shape[1]):
            blocksQ[i][j]= j%len(cols_B_mod)
    Q_ass=assemble_block((QM,QV,blocksQ))




    B_augm=aux.col_echelon(np.concatenate([B_ass, Q_ass]),field)[0]
    B_ech , Q_ech = B_augm[:B_ass.shape[0]],B_augm[B_ass.shape[0]:]

    blocks_A=[[ M[i][j]*VxPy_full[b] for j,b in enumerate(blocks_i)]for i,blocks_i in enumerate(blocks)]
    blocks_P= [[ PM[i][j]*PV[b] for j,b in enumerate(blocks_i)]for i,blocks_i in enumerate(blocksP)]
    for y,j in enumerate(cols_B):
        for i in range(blocks.shape[1],blocks.shape[0]):
            (blocks_A[i][j]) [:,zero_col_top_V[j%len(V)]] = M[i][j]*VxPy_red[blocksB[i-blocks.shape[1] ][y]]
        for i in range(blocks.shape[1]) :
            (blocks_P[i][j]) [:,zero_col_top_V[j%len(V)]] = QM[i][y]* QV[blocksQ[i ][y]]
    pass
     
    return np.concatenate([np.block(blocks_A),np.block(blocks_P)])
   



<<<<<<< HEAD
def row_echelon_block(X,V_ech,PV,field=Field.Field("Q")):
=======
def row_echelon_block(X,V_ech,PV):
    field=Field.Field("Q")
>>>>>>> main
    (M,V,blocks)=X
    A_block_ech= assemble_block((M,V_ech,blocks))
    diag_P= [ PV[y%len(V)] for y in range(blocks.shape[0])]
    P_block_ech=Field.build_block_diag_l(diag_P,field)
    return aux.row_echelon(np.concatenate([A_block_ech,P_block_ech],axis=1), field, max_col=A_block_ech.shape[1])