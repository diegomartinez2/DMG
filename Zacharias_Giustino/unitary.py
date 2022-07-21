import numpy as np
import math
import random
import cmath
from math import pi,e,log
y=60
z=60
UA=np.matrix([[1, 0, 0, 0], [0, math.cos(4*y), math.sin(4*y), 0], [0, math.sin(4*y), -(math.cos(4*y)), 0], [0, 0, 0, -1]])
UB=np.matrix([[1, 0, 0, 0], [0, math.cos(4*z), math.sin(4*z), 0], [0, math.sin(4*z), -(math.cos(4*z)), 0], [0, 0, 0, -1]])
IV=np.matrix('1 ;1 ;0 ;0')
X=np.matrix('1 ;0 ;1 ;0')
print (X)
one=UA*IV
print ("UAIV",one)
two=UB*one
print ("UBUAIV",two)
UAT=np.transpose(UA.real)
three=(UAT*two)
print("UAT(UBUA(IV))=UB(IV)",three)
a = np.squeeze(np.asarray(three))
b = np.squeeze(np.asarray(IV))
four=a*b
print ("UB(IV)(X)",four)
IVT=np.transpose(IV.real)
c=np.squeeze(np.asarray(four))
d=np.squeeze(np.asarray(IVT))
five=c*d
print ("IVT(UB(IV)(X))=UB(X)",five)
UBT=np.transpose(UB.real)
e=np.squeeze(np.asarray(UBT))
f=np.squeeze(np.asarray(five))
six=e*f
print ("UBTUB(X)=X",six)

#Reference: https://www.physicsforums.com/threads/unitary-transformation-using-python.951680/
