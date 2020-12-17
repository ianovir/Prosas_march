#============================================================
#Copyright (c) [2020] [Sebastiano Campisi - ianovir.com]
#
#Licensed to the Apache Software Foundation (ASF) under one
#or more contributor license agreements.  See the NOTICE file
#distributed with this work for additional information
#regarding copyright ownership.  The ASF licenses this file
#to you under the Apache License, Version 2.0 (the
#"License"); you may not use this file except in compliance
#with the License.  You may obtain a copy of the License at
#
#  http://www.apache.org/licenses/LICENSE-2.0
#
#Unless required by applicable law or agreed to in writing,
#software distributed under the License is distributed on an
#"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
#KIND, either express or implied.  See the License for the
#specific language governing permissions and limitations
#under the License.
#============================================================

from fractions import Fraction

def gcd(a, b):
    while b:
        a, b = b, a%b
    return a

def transposeMatrix(m):
    return map(list,zip(*m))

def getMatrixMinor(m,i,j):
    return [row[:j] + row[j+1:] for row in (m[:i]+m[i+1:])]

def getDeterminant(m):
    if len(m) == 2:
        return m[0][0]*m[1][1]-m[0][1]*m[1][0]
    det = Fraction(0)
    for c in range(len(m)):
        det += ((-1)**c)*m[0][c]*getDeterminant(getMatrixMinor(m,0,c))
    return det

def invertMatrix(m):
    det = getDeterminant(m)
    if len(m)==1: return m
    if len(m) == 2:
        return [[m[1][1]/det, -1*m[0][1]/det],
                [-1*m[1][0]/det, m[0][0]/det]]
    cof = []
    for r in range(len(m)):
        crow = []
        for c in range(len(m)):
            minDet = getDeterminant(getMatrixMinor(m,r,c))
            crow.append(((-1)**(r+c)) * minDet)
        cof.append(crow)
    cof = transposeMatrix(cof)
    for r in cof:
        for c in cof:
            if det!=0:
                cof[r][c] = cof[r][c]/det
            else:
                return None
    return cof

def getIdentity(size):
  I = []
  for r in range(size):
    ir = []
    for c in range(size):
      ir.append(Fraction(1) if c==r else Fraction(0))
    I.append(ir)
  return I

def subtractMatrices(a, b):
  res = [] # todo: check input's dimensions
  for r in range(len(a)):
    rr = []
    for c in range(len(a[r])):
      rr.append(a[r][c]-b[r][c])
    res.append(rr)
  return res
  
def multiplyMatrices(a, b):
    zip_b = list(zip(*b))
    return [[sum(ele_a*ele_b for ele_a, ele_b in zip(row_a, col_b)) 
             for col_b in zip_b] for row_a in a]

def getAbsorbingStates(m):
    """
    Returns a list with the indices of the absorbing states of the matrix m
    """
	
    a=[]
    for r in range(len(m)):
      if(sum(m[r])==0): a.append(r)
    return a

def getTransientStates(m):
    """Returns a list with the indices of the non-absorbing states of the matrix m"""
    t= list(range(len(m)))
    for r in range(len(m)):
      if(sum(m[r])==0): t.remove(r)
    return t

def normalizeProbabilityMatrix(m):
    n = []
    for r in range(len(m)):
      nrow = []    
      sumr = max(sum(m[r]), 1)       
      for c in range(len(m)):              
        nrow.append(Fraction(m[r][c], sumr))
      n.append(nrow)
    return n

def getQ(m, t):
    """
    Returns the matrix Q as in the canonical form [[Q,R],[0,I]]
    m is the input matrix, t is the list of the transient states
    """

    Q = []
    for r in range(len(t)):
      qrow = []
      for c in range(len(t)):
        qrow.append(m[t[r]][t[c]])
      Q.append(qrow)  
    return Q

def getR(m, t, a):
    """
    Returns the matrix R as in the canonical form [[Q,R],[0,I]]
    m is the input matrix, t is the list of the transient states,
    a is the list of the absorbing states.
    """

    R = []
    for r in range(len(t)):
      qrow = []
      for c in range(len(a)):
        qrow.append(m[t[r]][a[c]])
      R.append(qrow)  
    return R

def lcm(x, y):
    return x * y // gcd(x, y)

def fraction_gcd(x, y):
    N = gcd(x.numerator, y.numerator)
    D = lcm(x.denominator, y.denominator)
    return Fraction(N, D)

def getGCD(f):
    res = f[0]
    for c in f[1::]:
      res = fraction_gcd(res , c)
    return res

def solve(m):
    """
	Input m is a list of list as transition matrix of the equivalent 
    Markov Chain; the first state (first row) is considered as the 
    starting one.
    
    Returns an array containing the probability of each absorbing state,
    None if there is no abs. state. If there are K absorbing states, the
    returned array will contain K+1 integer numbers, where the first K
    items are the numerators of the abs. states' probabilities, and the
    last one is the denominator.
    """
	
    #with the assumption that at least one terminal state is given:
    if(len(m)==2 or len(m)==1): return [1,1]
  
    #Normalizing the in. matrix and identifying the trans./abs. states:
    m = normalizeProbabilityMatrix(m)
    t = getTransientStates(m)
    a = getAbsorbingStates(m)
	
    if len(a) >0:
        print( str(len(a)) + " absorbing state" + ("" if len(a)<=1 else "s" ))
    else:
        print("No absorbing state detected")
        return
    
    #Getting the matrices Q and R as in the canonical form:
    Q = getQ(m,t)
    R = getR(m,t,a)
    I = getIdentity(len(Q))
    I_Q = subtractMatrices(I, Q)
    
    #Getting the fundamental matrix
    N = invertMatrix(I_Q)
    F = multiplyMatrices(N,R)
    
    #packing the result with a common denominator:
    gcd = getGCD(F[0]).denominator
    res=[]
    sum = 0
    for r in F[0]:
      val = int(r.numerator*(gcd/r.denominator))
      sum+=val
      res.append(val)
    res.append(sum)  
    return res
