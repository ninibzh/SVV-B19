from math import *
import numpy as np
from control.matlab import *
import matplotlib.pyplot as plt

# Citation 550 - Linear simulation

# xcg = 0.25 * c

# Stationary flight condition

hp0    =    5862.791773674935 *0.3048	      # pressure altitude in the stationary flight condition [m]
V0     =    0.514444444 * 177.10160557043372         # true airspeed in the stationary flight condition [m/sec]
alpha0 =    5.047869299072158 *np.pi/180        # angle of attack in the stationary flight condition [rad]
th0    =    3.487237402960242*np.pi/180        # pitch angle in the stationary flight condition [rad]

# Aircraft mass
m      =      6208.879543162097    # mass [kg] #mass at 30 min

# aerodynamic properties
e      =   0.8          # Oswald factor [ ]
CD0    =   0.04         # Zero lift drag coefficient [ ]
CLa    =   5.084          # Slope of CL-alpha curve [ ]

# Longitudinal stability
Cma    = -0.5626            # longitudinal stabilty [ ]
Cmde   =  -1.1642          # elevator effectiveness [ ]

# Aircraft geometry

S      = 30.00	          # wing area [m^2]
Sh     = 0.2 * S         # stabiliser area [m^2]
Sh_S   = Sh / S	          # [ ]
lh     = 0.71 * 5.968    # tail length [m]
c      = 2.0569	          # mean aerodynamic cord [m]
lh_c   = lh / c	          # [ ]
b      = 15.911	          # wing span [m]
bh     = 5.791	          # stabilser span [m]
A      = b ** 2 / S      # wing aspect ratio [ ]
Ah     = bh ** 2 / Sh    # stabilser aspect ratio [ ]
Vh_V   = 100.11566510739169 	          # [ ]
ih     = -2 * pi / 180   # stabiliser angle of incidence [rad]

# Constant values concerning atmosphere and gravity

rho0   = 1.2250          # air density at sea level [kg/m^3] 
lambd = -0.0065         # temperature gradient in ISA [K/m]
Temp0  = 288.15          # temperature at sea level in ISA [K]
R      = 287.05          # specific gas constant [m^2/sec^2K]
g      = 9.81            # [m/sec^2] (gravity constant)

# air density [kg/m^3]  
rho    = rho0 * pow( ((1+(lambd * hp0 / Temp0))), (-((g / (lambd*R)) + 1)))
W      = m * g            # [N]       (aircraft weight)

# Constant values concerning aircraft inertia

muc    = m / (rho * S * c)
mub    = m / (rho * S * b)
KX2    = 0.019
KZ2    = 0.042
KXZ    = 0.002
KY2    = 1.25 * 1.114

# Aerodynamic constants

Cmac   = 0                      # Moment coefficient about the aerodynamic centre [ ]
CNwa   = CLa                    # Wing normal force slope [ ]
CNha   = 2 * pi * Ah / (Ah + 2) # Stabiliser normal force slope [ ]
depsda = 4 / (A + 2)            # Downwash gradient [ ]

# Lift and drag coefficient

CL = 2 * W / (rho * V0 ** 2 * S)              # Lift coefficient [ ]
CD = CD0 + (CLa * alpha0) ** 2 / (pi * A * e) # Drag coefficient [ ]

# Stabiblity derivatives

CX0    = W * sin(th0) / (0.5 * rho * V0 ** 2 * S)
CXu    = -0.02792
CXa    = +0.47966		# Positive! (has been erroneously negative since 1993) 
CXadot = +0.08330
CXq    = -0.28170
CXde   = -0.03728

CZ0    = -W * cos(th0) / (0.5 * rho * V0 ** 2 * S)
CZu    = -0.37616
CZa    = -5.74340
CZadot = -0.00350
CZq    = -5.66290
CZde   = -0.69612

Cmu    = +0.06990
Cmadot = +0.17800
Cmq    = -8.79415

CYb    = -0.7500
CYbdot =  0     
CYp    = -0.0304
CYr    = +0.8495
CYda   = -0.0400
CYdr   = +0.2300

Clb    = -0.10260
Clp    = -0.71085
Clr    = +0.23760
Clda   = -0.23088
Cldr   = +0.03440

Cnb    =  +0.1348
Cnbdot =   0     
Cnp    =  -0.0602
Cnr    =  -0.2061
Cnda   =  -0.0120
Cndr   =  -0.0939


#Symmetric

C1 = np.array([[-2*muc*c, 0, 0, 0],
      [0, (CZa-2*muc)*c/V0, 0, 0],
      [0, 0, -c/V0, 0],
      [0, Cmadot*c/V0, 0, -2*muc*KY2]])
C2 = np.array([[CXu*V0, CXa, CZ0, CXq*V0/c],
      [CZu*V0, CZa, -CX0, V0/c*(CZq+2*muc)],
      [0, 0, 0, V0/c],
      [Cmu*V0, Cma, 0, Cmq*V0/c]])
C3 = np.array([[CXde], [CZde], [0], [Cmde]])
C_extra=C2+C1

A = (-np.linalg.inv(C1)).dot(C2)
B = (-np.linalg.inv(C1)).dot(C3)
C = np.array([[1, 0, 0, 0],
              [0, 1, 0, 0],
              [0, 0, 1, 0],
              [0, 0, 0, 1]])
D = (np.linalg.inv(C_extra)).dot(C3)
print(np.linalg.eigvals(A))

sys=ss(A,B,C,D)


#A-Symmetric (_a)

C1_a = np.array([[(CYbdot-2*mub)*b/V0, 0, 0, 0],
               [0, -0.5*b/V0, 0, 0],
               [0, 0, -8*mub*KX2, 8*mub*KXZ],
               [Cnbdot*b/V0, 0, 8*mub*KXZ, -8*mub*KZ2]])
C2_a = np.array([[CYb, CL, CYp*2*V0/b, (CYr-4*mub)*2*V0/b],
               [0, 0, (1-2*V0/b), 0],
               [Clb, 0, Clp*2*V0/b, Clr*2*V0/b],
               [Cnb, 0, Cnp*2*V0/b, Cnr*2*V0/b]])
C3_a = np.array([[CYda, CYdr],
               [0, 0],
               [Clda, Cldr],
               [Cnda, Cndr]])

C_extra_a=C2_a+C1_a


A_a = (-np.linalg.inv(C1_a)).dot(C2_a)
B_a = (-np.linalg.inv(C1_a)).dot(C3_a)
C_a = np.array([[1, 0, 0, 0],
              [0, 1, 0, 0],
              [0, 0, 1, 0],
              [0, 0, 0, 1]])
D_a =(np.linalg.inv(C_extra_a)).dot(C3)
print(np.linalg.eigvals(A_a))

'''
#Response for symmetric flight
x0=np.matrix([[V0],[alpha0], [th0], [0.05]]) 
t=np.arange(0.0,1000.01,1) 

#Defining inputs
w1=1
w2=2
u1 = np.sin(t*w1) 
u2 = np.sin(t*w2)

y1=initial(sys,t,x0)  

speed=[]
#print(len(t))
#print(len(y1[0]))

for i in range(len(y1[0])):
    speed.append(y1[0][i][0])
    
y1,tdum,xdum=lsim(sys,U=u1,T=t) 
y2,tdum,xdum=lsim(sys,U=u2,T=t)
plt.plot(t,speed)
plt.show()
plt.plot(t,y1,t,y2)
plt.plot()
'''

#starting tas, alt, pitch, aoa, time (sec), mass
#171.30657236187153 5901.040529995785 8.478089398438776 5.914765769730046 3237.0 6266.500016181866 #Phugoid
#177.10160557043372 5862.791773674935 3.487237402960242 5.047869299072158 3717.0 6208.879543162097 #Dutch Roll
#175.2228241367164 5853.042471866842 5.565115261585568 5.812438264793599 3635.0 6218.7233594426825 #Short Period

#eigen values


#phugoid
'''
(array([-1.28404348+1.89207559j, -1.28404348-1.89207559j,
        0.00547665+0.13477211j,  0.00547665-0.13477211j]), array([[ 1.02712457e-04-1.72145646e-04j,  1.02712457e-04+1.72145646e-04j,
        -3.54206367e-04+9.23701801e-03j, -3.54206367e-04-9.23701801e-03j],
       [ 7.53371008e-01+0.00000000e+00j,  7.53371008e-01-0.00000000e+00j,
        -1.09614251e-03+5.16459670e-02j, -1.09614251e-03-5.16459670e-02j],
       [ 5.56023761e-01-3.51096837e-01j,  5.56023761e-01+3.51096837e-01j,
         9.98622070e-01+0.00000000e+00j,  9.98622070e-01-0.00000000e+00j],
       [-2.70508693e-05+8.18692133e-04j, -2.70508693e-05-8.18692133e-04j,
         2.97932012e-06+7.33166340e-05j,  2.97932012e-06-7.33166340e-05j]]))

#Dutch Roll
(array([-4.25400348+0.j        , -0.33354818+2.04343478j,
       -0.33354818-2.04343478j, -0.01082982+0.j        ]), array([[-6.66815998e-02+0.00000000e+00j,  8.19777524e-01+0.00000000e+00j,
         8.19777524e-01-0.00000000e+00j, -1.44902343e-02+0.00000000e+00j],
       [-9.97140558e-01+0.00000000e+00j, -1.19074454e-01-5.59916090e-01j,
        -1.19074454e-01+5.59916090e-01j, -9.99894673e-01+0.00000000e+00j],
       [-3.54362697e-02+0.00000000e+00j, -9.89002645e-03+4.72517248e-04j,
        -9.89002645e-03-4.72517248e-04j, -9.04626569e-05+0.00000000e+00j],
       [-2.92289877e-03+0.00000000e+00j,  9.48042096e-04-1.34581890e-02j,
         9.48042096e-04+1.34581890e-02j, -8.16976316e-04+0.00000000e+00j]]))


#short period
(array([-1.32306705+1.94243638j, -1.32306705-1.94243638j,
        0.00344641+0.13572562j,  0.00344641-0.13572562j]), array([[ 8.55676543e-05-1.75691889e-04j,  8.55676543e-05+1.75691889e-04j,
        -2.05559787e-04+8.81802751e-03j, -2.05559787e-04-8.81802751e-03j],
       [ 7.53094383e-01+0.00000000e+00j,  7.53094383e-01-0.00000000e+00j,
         1.06459391e-03+5.01779749e-02j,  1.06459391e-03-5.01779749e-02j],
       [ 5.56225379e-01-3.51370876e-01j,  5.56225379e-01+3.51370876e-01j,
         9.98700772e-01+0.00000000e+00j,  9.98700772e-01-0.00000000e+00j],
       [-2.78082416e-05+8.04611744e-04j, -2.78082416e-05-8.04611744e-04j,
         1.79213204e-06+7.05773353e-05j,  1.79213204e-06-7.05773353e-05j]]))

'''

