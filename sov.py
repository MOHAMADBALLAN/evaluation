# -*- coding: utf-8 -*-
"""
Created on Sat Oct 20 12:32:11 2018

@author: mm
"""
import re
from math import *
import sys
pre_test=sys.argv[1] # the protein secondry structure conformations sequnces file(output of the predection method)
idss=sys.argv[2]    # ids of the blind dataset
S_PRE=[]
S_OB=[]
VH=[]
VE=[]
V_=[]
f=open(pre_test)
l=[]
for line in f:
    line= line.rstrip()
    line=line.split()
    if  len(line) >0:
      s1= (line[2] )
      S_PRE.append(s1)


l=[]
ids=open(idss)
for line in ids:
   l.append( line.rstrip())
#print l
for i in l:
    #k=[]
    ss=open("/home/mm/Desktop/evaluation/dssp/"+i+".dssp")   # observed protein secondry structure conformations(output of dssp ) 
    for line in ss :
       if line[0]!=">":
        s2=( line.rstrip())
        S_OB.append(s2)
#print len(S_PRE[10])
#print len(S_OB[10])

for S2,S1 in zip(S_PRE,S_OB):
    obs=[]
    pre=[]
    lo=[]
    lp=[]
    l=[]
    for i in S1:
        if i=="H":
           l.append(i)
           
    totalH= len(l)
   # print totalH
    
    
    for match in re.finditer("H+",S1):
        lo.append( ( match.start(), match.end()))
     
        obs.append( match.group())

     #start = max(min1,min2)
   # end = min(max1,max2)
   # d = end - start
#print l

#print obs
    for match in re.finditer("H+",S2):
        pre.append( match.group())
        lp.append( ( match.start(), match.end()))
    sum_=0
    sum1=0
    for i in lo:
    
        for j in lp:
        
            start=max(min(i),min(j))
        #print start
        #print end
            end=min(max(i),max(j))
            minov=end -start
            if minov > 0:
           # print i,j, minov
           
               maxov=max(max(i),max(j))- min(min(i),min(j))
          
               segma = min( (maxov-minov),minov,((max(j)-min(j))/2),((max(i)-min(i))/2))
            
          
               sum_+= float (((minov)+(segma))/float(maxov))*float((max(i)-min(i)))
               sum1+=float((max(i)-min(i)))
    
    sumK=0
    for i ,j in zip(lo,lp):
         start=max(min(i),min(j))
         end=min(max(i),max(j))
         minov=end -start
         if minov < 0:
            sumK+=float(max(i)-min(i))
    NOH=sum1+sumK
   
       
                 
    #print sum_
          
           
            

            
    if totalH >0:
        
         SOVH=  (100/float(NOH))  * sum_
         VH.append ( SOVH)
    else:
        VH.append (0.0)
sum2_=0
for val in VH:
    sum2_+=val
F_SOVH =float(sum2_/len(VH))
print ("F_SOVH","=",F_SOVH)



#for E



for S2,S1 in zip(S_PRE,S_OB):
    obsE=[]
    prE=[]
    loE=[]
    lpE=[]
    lE=[]
    for i in S1:
        if i=="E":
           lE.append(i)
           
    totalE= len(lE)
   # print totalH
    
    
    for match in re.finditer("E+",S1):
        loE.append( ( match.start(), match.end()))
     
        obsE.append( match.group())

     #start = max(min1,min2)
   # end = min(max1,max2)
   # d = end - start
#print l

#print obs
    for match in re.finditer("E+",S2):
        prE.append( match.group())
        lpE.append( ( match.start(), match.end()))
    sumE_=0
    sum1=0
    for i in loE:
    
        for j in lpE:
        
            start=max(min(i),min(j))
        #print start
        #print end
            end=min(max(i),max(j))
            minov=end -start
            if minov > 0:
           # print i,j, minov
           
               maxov=max(max(i),max(j))- min(min(i),min(j))
          
               segma = min( (maxov-minov),minov,((max(j)-min(j))/2),((max(i)-min(i))/2))
            
          
               sumE_+= float (((minov)+(segma))/float(maxov))*float((max(i)-min(i)))
               sum1+=float(max(i)-min(i))
    sumK=0  
    for i ,j in zip(lo,lp):
        start=max(min(i),min(j))
        end=min(max(i),max(j))
        minov=end -start
        if minov < 0:
           sumK+=float(max(i)-min(i))
    NOE=sum1+sumK
   
          
           
            

            
    if totalE >0:
        
         SOVE=  (100/float(NOE))  * sumE_
         VE.append ( SOVE)
    else:
        VE.append (0.0)
sum2E_=0
for val in VE:
    sum2E_+=val
F_SOVE =float(sum2E_/len(VE))
print ("F_SOVE","=",F_SOVE)

for S2,S1 in zip(S_PRE,S_OB):
    obsc=[]
    prc=[]
    loc=[]
    lpc=[]
    lc=[]
    for i in S1:
        if i=="-":
           lc.append(i)
           
    totalc= len(lc)
   # print totalH
    
    
    for match in re.finditer("-+",S1):
        loc.append( ( match.start(), match.end()))
     
        obsc.append( match.group())

     #start = max(min1,min2)
   # end = min(max1,max2)
   # d = end - start
#print l

#print obs
    for match in re.finditer("-+",S2):
        prc.append( match.group())
        lpc.append( ( match.start(), match.end()))
    sumc_=0
    sum1=0
    for i in loc:
    
        for j in lpc:
        
            start=max(min(i),min(j))
        #print start
        #print end
            end=min(max(i),max(j))
            minov=end -start
            if minov > 0:
           # print i,j, minov
           
               maxov=max(max(i),max(j))- min(min(i),min(j))
          
               segma = min( (maxov-minov),minov,((max(j)-min(j))/2),((max(i)-min(i))/2))
            
          
               sumc_+= float (((minov)+(segma))/float(maxov))*float((max(i)-min(i)))
               sum1=float (max(i)-min(i))
    sumK=0  
    for i ,j in zip(lo,lp):
        start=max(min(i),min(j))
        end=min(max(i),max(j))
        minov=end -start
        if minov < 0:
           sumK+=float(max(i)-min(i))
    NO_=sum1+sumK
    #print sum_
          
           
            

            
    if totalc >0:
        
         SOVc=  (100/float(NO_))  * sumc_
         V_.append ( SOVc)
    else:
        V_.append (0.0)
sum2c_=0
for val in V_:
    sum2c_+=val
F_SOVC =float(sum2c_/len(V_))
print ("F_SOVC","=",F_SOVC)

    

               
    





      
       
        
        
      
      
      
    