# resmap
Calculate the secular resonant structure for a given configuration of
planets

Run using 

$> python sec_resonance_map.py input.file param.file

where input.file contains the initial conditions for the planetary
system. This file should contain the number of planets, their masses in
solar units, semi-major axes in AU, eccentricity, argument of
pericentre, inclination, and longitude of ascending node. An example of
appropriate input file for two planets is shown below

NUM 2
m 0.000954786 0.000285837
a 5.20245 9.554841
e 0.0474622 0.0575481
w 13.983865 88.719425
I 1.30667 2.48795
O 100.0381 113.1334 

This is included in the respository as "ic.dat".  The other file is the
parameter file for the code. An example, given in the repository as
"cntrl.dat" is shown here:

starting time=0.0
number of intervals=1
time interval=628.
radial zones=1000
inner radius=0.1
outer radius=30.
include GR=0
star mass=1.0
test e=0.0
test a=0.1
test I=0.0
test w=0.0
test O=0.0 

It only has the options shown here. "starting time" sets the zero point.
Just set to 0.0 unless you know what you are doing. The "number of
intervals" will usually be 1, but you can add many intervals to average
over.  "time interval" is used for "number of intervals" > 1. "radial
zones" gives the number of locations used for the test particle. "inner
radius" and "outer radius" is just that, in AU.  "include GR" 0 for no
and 1 for yes.  Not correctly implemented yet, so keep on 0. "star mass"
is in solar masses.  The "test" values give parameters chosen for the
test particle.  Only "e" and "I" really matter, as "a" will be set
according to the zones above.  These just serves to initialze other
values that may be used at some point. 

