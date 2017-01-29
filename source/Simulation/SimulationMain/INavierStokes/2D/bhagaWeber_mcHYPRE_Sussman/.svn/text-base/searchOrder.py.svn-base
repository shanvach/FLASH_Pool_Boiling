#!/usr/bin/env python
import re

logfileA=open("ApressureWrite.txt_001000")
logfileB=open("BpressureWrite.txt_001000")
logfileC=open("CpressureWrite.txt_001000")
logfileD=open("DpressureWrite.txt_001000")

a00 = []
a11 = []
a22 = []
a33 = []
a44 = []
a55 = []
a66 = []
a77 = []
a88 = []
b00 = []
b11 = []
b22 = []
b33 = []
b44 = []
b55 = []
b66 = []
b77 = []
b88 = []
c00 = []
c11 = []
c22 = []
c33 = []
c44 = []
c55 = []
c66 = []
c77 = []
c88 = []
d00 = []
d11 = []
d22 = []
d33 = []
d44 = []
d55 = []
d66 = []
d77 = []
d88 = []
alist=[]
blist=[]
clist=[]
dlist=[]

i=0
for line in logfileA:
	cols = [float(x) for x in line.split()]
	a0, a1, a2, a3, a4, a5, a6, a7, a8 = \
	   cols[0],cols[1], cols[2], cols[3], cols[4], cols[5], cols[6], cols[7], cols[8]
        a00.append(a0), a11.append(a1), a22.append(a2), a33.append(a3), a44.append(a4), \
       	   a55.append(a5), a66.append(a6), a77.append(a7), a88.append(a8)
        alist.append([a6, a7])
        #print i,alist[i][0],alist[i][1]
        i = i+1

i=0
for line in logfileB:
        cols = [float(x) for x in line.split()]
        b0, b1, b2, b3, b4, b5, b6, b7, b8 = \
           cols[0],cols[1], cols[2], cols[3], cols[4], cols[5], cols[6], cols[7], cols[8]
        b00.append(b0), b11.append(b1), b22.append(b2), b33.append(b3), b44.append(b4), \
           b55.append(b5), b66.append(b6), b77.append(b7), b88.append(b8)
        blist.append([b6, b7])
        #print i,blist[i][0],blist[i][1]
        i = i+1

i=0
for line in logfileC:
	cols = [float(x) for x in line.split()]
	c0, c1, c2, c3, c4, c5, c6, c7, c8 = \
           cols[0],cols[1], cols[2], cols[3], cols[4], cols[5], cols[6], cols[7], cols[8]
	c00.append(c0), c11.append(c1), c22.append(c2), c33.append(c3), c44.append(c4), \
           c55.append(c5), c66.append(c6), c77.append(c7), c88.append(c8)
        clist.append([c6, c7])
        #print i,c6,c7,c8
        i = i+1

i=0
for line in logfileD:
        cols = [float(x) for x in line.split()]
        d0, d1, d2, d3, d4, d5, d6, d7, d8 = \
           cols[0],cols[1], cols[2], cols[3], cols[4], cols[5], cols[6], cols[7], cols[8]
        d00.append(d0), d11.append(d1), d22.append(d2), d33.append(d3), d44.append(d4), \
           d55.append(d5), d66.append(d6), d77.append(d7), d88.append(d8)
        dlist.append([d6, d7])
        #print i,d6,d7,d8
        i = i+1

print '\n',"Lengths",len(a00), len(alist), alist[255][0],alist[255][1],'\n'

#+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_
#+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_

aSum = 0.0
aCount = 0.0

for j in range(0,len(alist)):
	for i in range(0,len(dlist)):
		if dlist[i][0] == alist[j][0]: 
			if dlist[i][1] == alist[j][1]: 
                		#print "HERE",j,i,"...",dlist[i][0],dlist[i][1],abs(a88[j]-d88[i])
				aSum   +=  (a88[j]-d88[i])**2.0
				aCount += 1

print '\n',"Finals for A",aSum,aCount,'\n'


#+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_
#+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_

bSum = 0.0
bCount = 0.0

for j in range(0,len(blist)):
        for i in range(0,len(dlist)):
                if dlist[i][0] == blist[j][0]:
                        if dlist[i][1] == blist[j][1]:
                                #print "HERE",j,i,"...",dlist[i][0],dlist[i][1],abs(b88[j]-d88[i])
                                bSum   +=  (b88[j]-d88[i])**2.0
                                bCount += 1

print '\n',"Finals for B",bSum,bCount,'\n'

#+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_
#+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_

cSum = 0.0
cCount = 0.0

for j in range(0,len(clist)):
        for i in range(0,len(dlist)):
                if dlist[i][0] == clist[j][0]:
                        if dlist[i][1] == clist[j][1]:
                                #print "HERE",j,i,"...",dlist[i][0],dlist[i][1],abs(c88[j]-d88[i])
                                cSum   +=  (c88[j]-d88[i])**2.0
                                cCount += 1

print '\n',"Finals for C",cSum,cCount,'\n'

#+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_
#+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_


print '\n',"Total AB",(1.0/aCount)*(aSum)**0.5,(1.0/bCount)*(bSum)**0.5
print '\n',"Total AB",(1.0/(aCount)**0.5)*(aSum)**0.5,(1.0/(bCount)**2.0)*(bSum)**0.5
print '\n',"Total BC",(1.0/bCount)*(bSum)**0.5,(1.0/cCount)*(cSum)**0.5
print '\n',"Total BC",(1.0/(bCount)**0.5)*(bSum)**0.5,(1.0/(cCount)**2.0)*(cSum)**0.5

