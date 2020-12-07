## the following python code find the 192 newkets using 192 group elements for each original kets and also find the corresponding phases for each newkets
import numpy as np

#load the sites permutation under 192 group elements

d = np.load("site transformation.npy")

# load the 16 sites' position

posits = np.loadtxt("16_positions.tsv")

b = posits.reshape(16,3)

## construct 192 16*16 matrix indicate site permutations

rep = np.zeros((192,16,16),dtype=int)

for i in range (0,192):

    for j in range (0,16):

        for k in range (0,16):

            if np.array_equal(d[i][j],b[k]) == True:

                rep[i][j][k]=1   

def bin_array(num, m):

#Convert a positive integer num into an m-bit bit vector

    return np.array(list(np.binary_repr(num).zfill(m))).astype(np.int8)

#define the original kets

ket = np.zeros((65536,16), dtype = int)

for i in range(0,65536):

    ket[i]=bin_array(i,16) 

# find out the rotation on local axis

#d3h = (id,c3z,c3z2,c2y,c2yc3z,c2yc3z2,inv,invc3z,invc3z2,invc2y,invc2yc3z,invc2yc3z2)

effect_op2 = np.loadtxt("effect_op2.tsv", dtype=int)

# define new_ket from the result of 192 operations act on the original kets




new_ket=np.zeros((65536,192,16), dtype = int)

for i in range (0,65536):

    for j in range (0,192):

        new_ket[i][j]=np.matmul(rep[j], ket[i])

# flip the spin under C2y

        if (effect_op2[j]==4 or effect_op2[j]==5 or effect_op2[j]==6 or effect_op2[j]==10 or effect_op2[j]==11 or effect_op2[j]==12):

            new_ket[i][j]=1-new_ket[i][j]

# convert binary to decimal

def binary_to_decimal(binary):

    decimal=binary[0]*2**15+binary[1]*2**14+binary[2]*2**13+binary[3]*2**12+binary[4]*2**11+binary[5]*2**10+binary[6]*2**9+binary[7]*2**8+binary[8]*2**7+binary[9]*2**6+binary[10]*2**5+binary[11]*2**4+binary[12]*2**3+binary[13]*2**2+binary[14]*2**1+binary[15]*2**0

    return decimal


#define newket, it store all non zero position on each new basis kets

newket=np.zeros((65536,192), dtype = int)

for i in range (0,65536):

    for j in range (0,192):

        newket[i][j]=binary_to_decimal(new_ket[i][j])

# define phase 
phase = np.zeros((65536,192),dtype=complex)

for i in range (0,65536):
    for j in range (0,192):
        
        if (effect_op2[j]==1 or effect_op2[j]==7):
            phase[i][j]=1

        elif (effect_op2[j]==2 or effect_op2[j]==8):
            ele=1    
            for k in range(0,16):
                ele=ele*np.exp(-1j*(new_ket[i][j][k]-1/2)*2*np.pi/3)
                phase[i][j]=ele

        elif (effect_op2[j]==3 or effect_op2[j]==9):
            ele=1
            for k in range(0,16):
                ele=ele*np.exp(+1j*(new_ket[i][j][k]-1/2)*2*np.pi/3)
                phase[i][j]=ele 

        elif (effect_op2[j]==4 or effect_op2[j]==10):

        ###  c2y will give 1 or -1, count the number if even is 1, if odd is -1.
            ele=1    
            for k in range(0,16):
                ele=ele*(2*(new_ket[i][j][k]-1/2))
                phase[i][j]=ele

        elif (effect_op2[j]==5 or effect_op2[j]==11):
            ele1=1
            ele2=1
            for k in range(0,16):
                ele1=ele1*(2*(new_ket[i][j][k]-1/2))
                ele2=ele2*np.exp(-1j*(new_ket[i][j][k]-1/2)*2*np.pi/3)
                phase[i][j]=ele1*ele2

        elif (effect_op2[j]==6 or effect_op2[j]==12):
            ele1=1
            ele2=1
            for k in range(0,16):
                ele1=ele1*(2*(new_ket[i][j][k]-1/2))
                ele2=ele2*np.exp(+1j*(new_ket[i][j][k]-1/2)*2*np.pi/3)
                phase[i][j]=ele1*ele2

np.save('newket',newket)
np.save('phase',phase)

