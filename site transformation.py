import numpy as np
#First load the 48 point group elements, stored as a 4X3 matrices.
#the first 3X3 components are a rotation, then a translation is given by a 3-component vector 
#in the fourth column
symfuns = np.loadtxt("symm_elem.3.tsv")
symfuns = symfuns.reshape(48,3,4)
a = np.zeros((192,3,4))


for count in range(0,48):
  a[count] = symfuns[count]
  a[count + 1*48] = symfuns[count]
  a[count + 2*48] = symfuns[count]
  a[count + 3*48] = symfuns[count]
#2nd make 7 copies of the point group elements, with extra translations
#unlike bilbao, we assume rotation is performed first, then the translation.


for count in range(0,48):
   a[48 + count][1][3] += 0.5
   a[48 + count][2][3] += 0.5
   a[2*48 + count][0][3] += 0.5
   a[2*48 + count][2][3] += 0.5
   a[3*48 + count][0][3] += 0.5
   a[3*48 + count][1][3] += 0.5


##16 sites shift_adj
def shift_adj(opA):
   opc=opA
   if opc[2][3] >= 1:        
        opc[2][3] += -1
   if opc[2][3] < 0:
        opc[2][3] += 1
   if opc[1][3] >= 1:
        opc[1][3] += -1
   if opc[1][3] < 0:
        opc[1][3] += 1
   if opc[0][3] >= 1:
      opc[0][3] += -1 
   if opc[0][3] < 0:
      opc[0][3] += 1
   return opc


for count in range(0,48*4):
   a[count]= shift_adj(a[count])

ele= np.zeros((192*3,4),dtype=int)
ele=a.reshape(192*3,-1)


posits = np.loadtxt("16_positions.tsv")

b = posits.reshape(16,3)
c = np.zeros((192,16),dtype=int)

# 16 sites function
def vec_shift_adj(opA):
   opc=opA
   if opc[2] >= 1:
          opc[2] += -1
   if opc[2] < 0:
          opc[2] += 1
   if opc[1] >= 1:
        opc[1] += -1
   if opc[1] < 0:
        opc[1] += 1
   if opc[0] >= 1:
        opc[0] += -1
   if opc[0] < 0:
        opc[0] += 1
   return opc

d = np.zeros((192,16,3))

for count5 in range (0,192):
   for count1 in range (0,16):
      for count2 in range (0,3):
         for count3 in range (0,3):
            d[count5][count1][count2] += a[count5][count2][count3]*b[count1][count3]             
         d[count5][count1][count2] += a[count5][count2][3]
      d[count5][count1]=vec_shift_adj(d[count5][count1])


for count5 in range (0,192):
   for count1 in range (0,16):
      for count4 in range (0,16):
         if np.array_equal(d[count5][count1],b[count4]) == True:
            c[count5][count1] = count4+1
            
np.save('site transformation',d)
