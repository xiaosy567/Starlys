
#!/usr/bin/env python
'''
rotate abaout an arbitrary axis (Line)
'''

import numpy as np
import math
import sys
import os
from math import cos, sin, sqrt, acos
from numpy import linalg, dot
from numpy import array
from numpy import matrix

def rotate(axis_a1, axis_a2, point_martix, theta):
    '''
    This function calculate the coordinate after rotation about an 
    arbitrary axis (line).

    input: axis point 1 (start point);
           axis point 2 (end point);
           point martix np.array([ [x1, y1, z1, 1.0], [x2, y2, z2, 1.0], [x3, y3, z3, 1.0], ... ]);
           theta (rotation angle);
    '''
    # print "@rotate about an arbitrary axis line"
    # print axis_a1
    # print axis_a2
    # print point_martix
    # print theta

    #theta_rad=theta*3.1415926/180.0
    theta_rad=theta
    A=axis_a2[0]-axis_a1[0]
    B=axis_a2[1]-axis_a1[1]
    C=axis_a2[2]-axis_a1[2]
    L=sqrt(A*A + B*B + C*C)
    V=sqrt(B*B + C*C)
    # print A, B, C, L, V
    Do = np.array([ [1, 0, 0, -axis_a1[0]], 
                    [0, 1, 0, -axis_a1[1]], 
                    [0, 0, 1, -axis_a1[2]], 
                    [0, 0, 0, 1] ])
    Rx = np.array([ [1, 0, 0, 0],
                    [0, C/V, -B/V, 0],
                    [0, B/V, C/V, 0],
                    [0, 0, 0, 1] ])
    Ry = np.array([ [V/L, 0, -A/L, 0],
                    [0, 1, 0, 0],
                    [A/L, 0, V/L, 0],
                    [0, 0, 0, 1]  ])
    Rz = np.array([ [cos(theta_rad), -sin(theta_rad), 0, 0],
                    [sin(theta_rad), cos(theta_rad), 0, 0],
                    [0, 0, 1, 0],
                    [0, 0, 0, 1]  ])
    Ryr = np.array([ [V/L, 0, A/L, 0],
                    [0, 1, 0, 0],
                    [-A/L, 0, V/L, 0],
                    [0, 0, 0, 1]  ])
    Rxr = np.array([ [1, 0, 0, 0],
                    [0, C/V, B/V, 0],
                    [0, -B/V, C/V, 0],
                    [0, 0, 0, 1] ])
    Dor = np.array([ [1, 0, 0, axis_a1[0]], 
                    [0, 1, 0, axis_a1[1]], 
                    [0, 0, 1, axis_a1[2]], 
                    [0, 0, 0, 1] ])

    #print Do, '\n', matrix(point_martix).T
    P_Do=matrix(Do)*(matrix(point_martix).T)
    P_Rx=matrix(Rx)*(matrix(P_Do))
    P_Ry=matrix(Ry)*(matrix(P_Rx))
    P_Rz=matrix(Rz)*(matrix(P_Ry))
    P_Ryr=matrix(Ryr)*(matrix(P_Rz))
    P_Rxr=matrix(Rxr)*(matrix(P_Ryr))
    P_Dor=matrix(Dor)*(matrix(P_Rxr))
    # print P_Dor
    P_t=np.array(matrix(P_Dor).T)
    return P_t[0]
    
    #########################
    # The following methodology derive the same result as above.
    #########################
    #Trans_matrix=matrix(Dor)*matrix(Rxr)*matrix(Ryr)*matrix(Rz)*matrix(Ry)*matrix(Rx)*matrix(Do)
    #P2 = matrix(Trans_matrix)*matrix(point_martix).T
    #print np.array(Trans_matrix)
    #P_t=np.array(matrix(P2).T)
    #print P_t
    #return P_t

def angle_vectors(v1, v2):
    '''
    calculate the angle of two vectors.
    input: v1: vector-1; np.array([ x, y, z ])
           v2: vector-2;
           point martix np.array([ [x1, y1, z1, 1.0], [x2, y2, z2, 1.0], [x3, y3, z3, 1.0], ... ]);
           theta (rotation angle);
    '''
    pv1=matrix(v1)*(matrix(v1).T)
    v1u=np.array(v1)/(1.0*sqrt(pv1[0]))

    pv2=matrix(v2)*(matrix(v2).T)
    v2u=np.array(v2)/(1.0*sqrt(pv2[0]))

    p=matrix(v1u)*(matrix(v2u).T)

    theta=math.acos(p)
    # print theta, theta*180/3.1415926
    return theta

def return_unitvector(v1):
    pv1=matrix(v1)*(matrix(v1).T)
    v1u=np.array(v1)/(1.0*sqrt(pv1[0]))
    return v1u

# if __name__=="__main__":
   # refa1=np.array([6, -2, 0])
   # refa2=np.array([12, 8, 0])
   # p=np.array([ [3, 5, 0, 1], [10, 6, 0, 1], [1, 1, 0, 1], [3, 5, 0, 1], [5, 5, 2, 1] ])
   # theta=60.0
   # pt=rotate(refa1, refa2, p, theta)
   
   # refa1=np.array([1, 0, 0])
   # refa2=np.array([0, 1, 0])
   # print angle_vectors(refa1, refa2)

