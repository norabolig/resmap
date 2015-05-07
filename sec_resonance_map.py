import numpy as np
import math
import sys
from numpy import linalg

## MASS OF THE STAR ###
# Solar units
MC=0.0
#######################

## USED FOR INTEGRATION CONTROL IN CALC OF BINT ##
NSIMP=100
##################################################

## NOT EVERTHING BELOW IS USED.  EXPANSION FOR LATER ##

def load_ic(arr,num):
  ncheck=len(arr)
  if not ncheck==num+1: 
     print("SERIOUS ERROR IN ICS. ASSERTING.")
     sys.exit()
  tmp_arr=np.zeros(num)
  for ic in xrange(1,num+1): tmp_arr[ic-1]=float(arr[ic])
  return tmp_arr

def read_input(filename):
  fhndle=open(filename,"r")
  num=0
  for line in fhndle:
     if line[0:3]=="NUM":
         junk,cnum=line.split()
         num=int(cnum)
         m=np.zeros(num)
         a=np.zeros(num)
         e=np.zeros(num)
         w=np.zeros(num)
         I=np.zeros(num)
         O=np.zeros(num)
     if num==0:continue
     item=line.split()
     if   item[0]=="m": m=load_ic(item,num)
     elif item[0]=="a": a=load_ic(item,num)
     elif item[0]=="e": e=load_ic(item,num)
     elif item[0]=="w": w=load_ic(item,num)
     elif item[0]=="I": I=load_ic(item,num)
     elif item[0]=="O": O=load_ic(item,num)
     else: continue

  fhndle.close()

  return m,a,e,w,I,O 

def read_controls(filename):
        dict={"starting time":0.0,
              "number of intervals":1,
              "time interval":100*math.pi*2.0,
              "radial zones":1000,
              "inner radius":0.1,
              "outer radius":30.,
              "include GR":0, 
              "star mass":1.0,
              "test e":0.0,
              "test a":0.1,
              "test I":0.0,
              "test w":0.0,
              "test O":0.0}

        inlist=[]
        fhndle=open(filename,"r")
        for invalue in fhndle:
                k,v=invalue.rstrip().split("=")
                inlist.append(k)
                inlist.append(v)
        fhndle.close()

        for item in xrange(len(inlist)):
                if inlist[item]  =="starting time"       : dict[inlist[item]]=float(inlist[item+1])
                elif inlist[item]=="number of intervals" : dict[inlist[item]]=int  (inlist[item+1])
                elif inlist[item]=="time interval"       : dict[inlist[item]]=float(inlist[item+1])
                elif inlist[item]=="radial zones"        : dict[inlist[item]]=int  (inlist[item+1])
                elif inlist[item]=="inner radius"        : dict[inlist[item]]=float(inlist[item+1])
                elif inlist[item]=="outer radius"        : dict[inlist[item]]=float(inlist[item+1])
                elif inlist[item]=="include GR"          : dict[inlist[item]]=int  (inlist[item+1])
                elif inlist[item]=="star mass"           : dict[inlist[item]]=float(inlist[item+1])
                elif inlist[item]=="test e"              : dict[inlist[item]]=float(inlist[item+1])
                elif inlist[item]=="test a"              : dict[inlist[item]]=float(inlist[item+1])
                elif inlist[item]=="test I"              : dict[inlist[item]]=float(inlist[item+1])
                elif inlist[item]=="test w"              : dict[inlist[item]]=float(inlist[item+1])
                elif inlist[item]=="test O"              : dict[inlist[item]]=float(inlist[item+1])

        return dict


def nmm(a,m):
  return math.sqrt( (m+MC)/a**3 )

def simpson(y0,y1,y2,h):
  a=y0
  b=(3.*(y2-a)-4.*(y2-y1))/h
  c=((y2-a)/h-b)/h
  return h*(a+h*(0.5*b+h*c/3.))

def laplace_integrand(j,s,a,phi):
  return math.cos(j*phi)/((1.-2.*a*math.cos(phi)+a*a)**s)

def bint(j,s,a):
  dphi=2.0*math.pi/float(NSIMP)
  val=0.0
  for i in xrange(NSIMP):
     phi0=i*dphi
     phi1=phi0+0.5*dphi
     phi2=phi0+dphi
     y0=laplace_integrand(j,s,a,phi0)
     y1=laplace_integrand(j,s,a,phi1)
     y2=laplace_integrand(j,s,a,phi2)
     val+=simpson(y0,y1,y2,dphi)
  return val/math.pi

def ajk(a,b):
  if a>b: return b/a
  else: return a/b

def ajkbar(a,b):
  if a>b: return 1.0
  else: return a/b

def AijBij(m,a,n):
  N=len(m)
  A=np.zeros((N,N))
  B=np.zeros((N,N))
  # yeah, the following ain't pretty.  I'll fix this when speed is needed.
  for i in xrange(N):
    for k in xrange(i):
      factor=0.25*n[i]*(m[k]/(MC+m[i]))*ajk(a[i],a[k])*ajkbar(a[i],a[k])
      b1=bint(1.,1.5,ajk(a[i],a[k]))
      b2=bint(2.,1.5,ajk(a[i],a[k]))
      A[i][i]+=factor*b1
      A[i][k]=-factor*b2
#      A[i][i]+=0.25*n[i]*(m[k]/(MC+m[i]))*ajk(a[i],a[k])*ajkbar(a[i],a[k])*bint(1.,1.5,ajk(a[i],a[k]))
#      A[i][k]=-0.25*n[i]*(m[k]/(MC+m[i]))*ajk(a[i],a[k])*ajkbar(a[i],a[k])*bint(2.,1.5,ajk(a[i],a[k]))

      B[i][i]-=factor*b1
      B[i][k]=+factor*b1
#      B[i][i]-=0.25*n[i]*(m[k]/(MC+m[i]))*ajk(a[i],a[k])*ajkbar(a[i],a[k])*bint(1.,1.5,ajk(a[i],a[k]))
#      B[i][k]=+0.25*n[i]*(m[k]/(MC+m[i]))*ajk(a[i],a[k])*ajkbar(a[i],a[k])*bint(1.,1.5,ajk(a[i],a[k]))
    for k in xrange(i+1,N):
      factor=0.25*n[i]*(m[k]/(MC+m[i]))*ajk(a[i],a[k])*ajkbar(a[i],a[k])
      b1=bint(1.,1.5,ajk(a[i],a[k]))
      b2=bint(2.,1.5,ajk(a[i],a[k]))
      A[i][i]+=factor*b1
      A[i][k]=-factor*b2
      B[i][i]-=factor*b1
      B[i][k]=+factor*b1
  return A,B

def AB_part(m,a,apart,npart):
  N=len(m)
  A=0.
  B=0.
  for k in xrange(N):
      A+=0.25*npart*m[k]/(MC)*ajk(apart,a[k])*ajkbar(apart,a[k])*bint(1.,1.5,ajk(apart,a[k]))
  B=-1.0*A
  return A,B

def AjBj_part(m,a,apart,npart):
  N=len(m)
  A=np.zeros(N)
  B=np.zeros(N)
  for k in xrange(N):
      factor=0.25*npart*m[k]/(MC)*ajk(apart,a[k])*ajkbar(apart,a[k])
      A[k]=-factor*bint(2.,1.5,ajk(apart,a[k]))
      B[k]=+factor*bint(1.,1.5,ajk(apart,a[k]))

  return A,B

def dRdk(i,A,h,k):
   N=len(h)
   val=A[i][i]*k[i]
   for l in xrange(i):     val+=A[i][l]*k[l]
   for l in xrange(i+1,N): val+=A[i][l]*k[l]
   return val

def dRdh(i,A,h,k):
   N=len(h)
   val=A[i][i]*h[i]
   for l in xrange(i):     val+=A[i][l]*h[l]
   for l in xrange(i+1,N): val+=A[i][l]*h[l]
   return val

def dRdp(i,B,p,q):
   N=len(p)
   val=B[i][i]*p[i]
   for l in xrange(i):     val+=B[i][l]*p[l]
   for l in xrange(i+1,N): val+=B[i][l]*p[l]
   return val

def dRdq(i,B,p,q):
   N=len(q)
   val=B[i][i]*q[i]
   for l in xrange(i):     val+=B[i][l]*q[l]
   for l in xrange(i+1,N): val+=B[i][l]*q[l]
   return val


def hdot(i,A,h,k):
   return dRdk(i,A,h,k)

def kdot(i,A,h,k):
   return -dRdh(i,A,h,k)

def pdot(i,B,p,q):
   return dRdq(i,B,p,q)
 
def qdot(i,B,p,q):
   return -dRdp(i,B,p,q)

def set_hk(e,w):
   N=len(e)
   h=np.zeros(N)
   k=np.zeros(N)
   for i in xrange(N): h[i]=e[i]*math.sin(w[i])
   for i in xrange(N): k[i]=e[i]*math.cos(w[i])
   return h,k

def set_pq(I,O):
   N=len(I)
   p=np.zeros(N)
   q=np.zeros(N)
   for i in xrange(N): p[i]=I[i]*math.sin(O[i])
   for i in xrange(N): q[i]=I[i]*math.cos(O[i])
   return p,q

def get_ew(h,k):
   N=len(h)
   e=np.zeros(N)
   w=np.zeros(N)
   for i in xrange(N):
     w[i]=math.atan2(h[i],k[i]) 
     if w[i]<0.0:w[i]+=math.pi*2.0
     sinw=math.sin(w[i])
     if sinw>0.1 or sinw < -0.1:e[i]=h[i]/sinw
     else:
       e[i]=k[i]/math.cos(w[i])
   return e,w

def get_IO(p,q):
   N=len(p)
   I=np.zeros(N)
   O=np.zeros(N)
   for i in xrange(N):
     O[i]=math.atan2(p[i],q[i]) 
     if O[i]<0.0:O[i]+=math.pi*2.0
     sinO=math.sin(O[i])
     if sinO>0.1 or sinO < -0.1:I[i]=p[i]/sinO
     else:
       I[i]=q[i]/math.cos(O[i])
   return I,O

def hkdots_GR(h,k,a,n):
   N=len(h)
   hdotgr=np.zeros(N)
   kdotgr=np.zeros(N)
   c=10070.0
   wgrfac=3.*MC/c**2
   for i in xrange(N):
     wgr=wgrfac*n[i]/a[i]
     hdotgr[i]= k[i]*wgr
     kdotgr[i]=-h[i]*wgr
   return hdotgr,kdotgr

def main():

   if not len(sys.argv)==3:
     print("\nUSAGE:")
     print("python ./sec_resonance_map.py input_file param_file\n")
     print("Please give input file for initial conditions and parameter file for script control.\n")
     sys.exit()

   mIC,aIC,eIC,wIC,IIC,OIC=read_input(sys.argv[1])

   global MC

   dict=read_controls(sys.argv[2])
   T0=dict["starting time"]
   NTIME=dict["number of intervals"]
   DTIME=dict["time interval"]
   NR=dict["radial zones"]
   RIN=dict["inner radius"]
   ROUT=dict["outer radius"]
   INCLUDEGR=dict["include GR"]
   MC=dict["star mass"]
   etest_part=dict["test e"] 
   atest_part=dict["test a"] 
   itest_part=dict["test I"] 
   wtest_part=dict["test w"] 
   otest_part=dict["test O"] 
   mtest_part=0.0

### INITIALIZE ARRAYS ##
   N=len(mIC)
   m=np.zeros(N)
   a=np.zeros(N)
   e=np.zeros(N)
   I=np.zeros(N)
   w=np.zeros(N)
   O=np.zeros(N)
   n=np.zeros(N)
   twopi=2.0*math.pi
#######################

## RADIANS ####################
   for i in xrange(N):
     m[i]=mIC[i]
     a[i]=aIC[i]
     e[i]=eIC[i]
     I[i]=IIC[i]*math.pi/180.
     O[i]=OIC[i]*math.pi/180.
     w[i]=wIC[i]*math.pi/180.
     n[i]=nmm(a[i],m[i])
############################### 

   h,k=set_hk(e,w)
   p,q=set_pq(I,O) 

### As and Bs ###############################
# Get A and B MATRICES AND EIGENVALUES
   A,B=AijBij(m,a,n)
   lA,vA = linalg.eig(A)
   lB,vB = linalg.eig(B)
#############################################

## CALCULATE SECULAR THEORY!!! ##############

   #print "#A  ",A*2.*180.
   #print "#B  ",B*2.*180.
   #print "#lA ",lA*2.*180.
   #print "#lB ",lB*2.*180.
   #print "vA ",vA
   #print "vB ",vB
   #print "h ",h
   #print "k ",k
   #print "p ",p
   #print "q ",q

   vAinv = linalg.inv(vA)
   vBinv = linalg.inv(vB)

   ssinb = vAinv.dot(h)
   scosb = vAinv.dot(k)

   tsing = vBinv.dot(p)
   tcosg = vBinv.dot(q)

   #print "SSINB ",ssinb
   #print "SCOSB ",scosb

   eBeta=np.zeros(N)
   eScale=np.zeros(N)
   iBeta=np.zeros(N)
   iScale=np.zeros(N)
     
     
   for i in xrange(N):
      eBeta[i]=math.atan2(ssinb[i],scosb[i])
      iBeta[i]=math.atan2(tsing[i],tcosg[i]) 
      sb=math.sin(eBeta[i])
      if sb == 0.0: eScale[i]=scosb[i]/math.cos(eBeta[i])
      else: eScale[i]=ssinb[i]/sb
      sb=math.sin(iBeta[i])
      if sb == 0.0: iScale[i]=tcosg[i]/math.cos(iBeta[i])
      else: iScale[i]=tsing[i]/sb
       
   #print "#eScale,eBeta ",eScale,eBeta*180./math.pi
   #print "#iScale,iBeta ",iScale,iBeta*180./math.pi

   ee = np.zeros((N,N))
   ii = np.zeros((N,N))

   for i in xrange(N):
    for j in xrange(N):
       ee[i][j]=eScale[j]*vA[i][j]
       ii[i][j]=iScale[j]*vB[i][j]
       #print " ee ",i,j,ee[i][j]

## NOW WE HAVE THE SECULAR THEORY. INTRODUCE TEST PARTICLE AT RANGE OF RADII ###

## RADIAL RANGE ##############################
## 
## PRINT HEADER
   print "# FOR TEST PARTICLE"
   print "# index, radius (AU), A (deg/yr), B (deg/yr), eforced, Iforced (deg), A eigenvalues (deg/yr) and B eigenvalues (deg/yr)"
##
## 
   for ir in xrange(NR):
     rad = RIN + (ROUT-RIN)*float(ir)/float(NR)

     for i in [N-1]:
       apart=rad*1.0
       npart=nmm(apart,mtest_part)
##############################################

     Apart,Bpart=AB_part(m,a,apart,npart)
     Aj,Bj=AjBj_part(m,a,apart,npart)
     # TEST PARTICLE A,B,Aj,Bj

     nu = np.zeros(N)
     mu = np.zeros(N)
     for i in xrange(N):
       for j in xrange(N):
         nu[i]+=Aj[j]*ee[j][i]
         mu[i]+=Bj[j]*ii[j][i]

     eforced=0.
     iforced=0.
     for itime in xrange(NTIME):
       time = T0 +itime*DTIME
       h0=0.0
       k0=0.0
       p0=0.0
       q0=0.0
       for i in xrange(N):
         h0+= -nu[i]/(Apart-lA[i])*math.sin(lA[i]*time+eBeta[i])
         k0+= -nu[i]/(Apart-lA[i])*math.cos(lA[i]*time+eBeta[i])
         p0+= -mu[i]/(Bpart-lB[i])*math.sin(lB[i]*time+iBeta[i])
         q0+= -mu[i]/(Bpart-lB[i])*math.cos(lB[i]*time+iBeta[i])
       eforced+=math.sqrt(h0**2+k0**2)
       iforced+=math.sqrt(p0**2+q0**2)

     eforced/=float(NTIME)
     iforced/=float(NTIME)



### OUTPUT!!!!
     #print "# FOR TEST PARTICLE"
     #print "# index, radius (AU), A (deg/yr), B (deg/yr), eforced, Iforced (deg), A eigenvalues (deg/yr) and B eigenvalues (deg/yr)"
     print ir,rad,Apart*2.*180.,Bpart*2.*180.,eforced,iforced,
     for i in xrange(N):print lA[i]*2.*180.,
     for i in xrange(N):print lB[i]*2.*180.,
     print ""

#############################################################################################################


main() 
