
import numpy as np
import matplotlib.pyplot as plt
from cfractions import Fraction
import scipy as sc
import networkx as nx
#maximum computation error
epsilon=1e-10

## COMPUTATIONS IN ALL FIELDS

#field F_2
class F_2():
    #definition
    is_one=False
    def __init__(self, input):
        try:#input is of type F_2
            self.is_one= input.is_one
        except:#input is of type bool/int
            self.is_one=bool(input%2)
    #override operators
    def __add__(self,b):
        return F_2(self.is_one^b.is_one)
    def __sub__(self,b):
        return self+b
    def __mul__(self,b):
        return F_2(self.is_one&b.is_one)
    def __truediv__(self,b):
        if not(b.is_one):
            raise ValueError("Divide by 0")
        else:
            return self
    def __neg__(self):
        return self
    def __eq__(self,b):
        try:
            return self.is_one==F_2(b).is_one
        except:
            raise ValueError("different types f")
    
    #display
    def __str__(self):
        return "c("+str(int(self.is_one))+")"
    def __repr__(self):
        return "c("+str(int(self.is_one))+")"

#display fractions as strings
#Fraction.__repr__=Fraction.__str__

predefined_to_Field= { "C":  lambda x: np.complex64(x) , "R":  lambda x: np.float64(x) , "Q":  lambda x: Fraction(x) , "F_2":  lambda x: F_2(x)  } 

#general field definition
class Field:
    descr='O'
    zero=0.
    one=1.
    to_Field=None
    def __init__(self,*args):
        # usual fields R, C, F_2, Q
        if len(args)==1:
            if args[0] in predefined_to_Field.keys():
                self.descr=args[0]
                self.to_Field= np.frompyfunc( predefined_to_Field[self.descr],1,1)
                self.zero= self.to_Field(0)
                self.one = self.to_Field(1)
            else:
                raise ValueError("No predefined field "+args[0])
        #user-defined fields given by 2 instances 0 and 1 of an object supporting +,*,-,:,==,... and a description 
        else:
            self.zero=args[0]
            self.one=args[1]
            self.descr=args[2]
            def to_Field(x):
                x=int(x)
                if x ==0:
                    return self.zero
                if x==1 :
                    return self.one
                if x<0:
                    return - to_Field(-x)
                return to_Field(x-1)

#Replacing numpy operations if the field is not R or C
#np.zeros


def bool_to_field(b,field):
    if b:
        return field.one
    return field.zero
       

def build_block_diag(M,N,field):
    return np.block([[M, field.to_Field(np.zeros((M.shape[0], N.shape[1]),dtype="i"))],[ field.to_Field(np.zeros((N.shape[0], M.shape[1]),dtype="i")),N]])



def build_block_diag_l(l,field):
    dim0= sum([ M.shape[0] for M in l])
    dim1= sum([ M.shape[1] for M in l])
    res= field.to_Field(np.zeros((dim0,dim1),dtype="i"))
    xi,yi=0,0
    for i,M in enumerate(l):
        res[xi:xi+M.shape[0],yi:yi+M.shape[1]] = M
        xi,yi=xi+M.shape[0],yi+M.shape[1]
    return res

            
def is_all_zero_mat(M,field):
    if len(M.flatten())==0:
        return True
    if field.descr in ['R','C']:
        return np.max(np.abs(M))<epsilon
    return all([m==field.zero for m in M.flatten()])

def is_all_zero_elem(x,field):
    return is_all_zero_mat(np.array([x]).reshape((1,1)),field)


#Remove computation errors
def flatten_zero(U,field):
    if field.descr in ['R','C']:
        V=U.flatten()
        V[np.where(np.abs(V)<epsilon)[0]]=0.
        return V.reshape(U.shape)
    return U




