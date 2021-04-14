# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 18:43:42 2018

@author: sana
"""

'''Generating input and data file to to the system---hairpin and tail attached to it---along with pore 
hairpin.xyz contains coordinates and atom type of the hairpin which is obtained from model-hairpin-only code'''



import numpy as np








'''conversion from real to lj unit'''

sigma = 2.5         #bead size in angstrom unit
epsilon = 0.2       # in kcal/mol

epsilonH = 1.0/epsilon    # strength of hydrophobic interaction        

massReal = 96.0
massReduced = 96.0/massReal                       # mass is 1.0 is reduced unit. It is taken each bead is of same mass i.e. 96 Da i.e. 96 gram/mole

#Bond coefficients

l0 = 3.0/sigma             # equilibrium bond length 4.6 angstrom unit
#kb = round((1./2.)*100.*epsilonH/l0**2.,3)     # 1/2 100 epsilonH /a^2
kb = 171.0 * (sigma**2)/epsilon   # harmonic constant for the harmonic potential. value is 171 kcal/mol/Angstrom^2 ( energy/distance^2)kb

#Angle coefficients

#kTheta = (1./2.)*20.*epsilonH      # in unit (energy/radian^2)
kTheta = 50./epsilon
theta0 = 105.                       # in degrees


#Dihedral coefficients 

#diheA = 1.0*epsilonH
#diheB = 1.6*epsilonH
#diheC = 1.0*epsilonH


#dielectric constant
dielectric = 80.0

#pair potential coordinates

sigmalj = 3.0
epsilonlj = 0.2
rlj = 2.5                        # global cutoff for the lj potential is taken as 2.5*sigma
rlocal= ((2.0)**(1.0/6.0))*sigmalj/sigma               # local cutoff for lj interaction to include only repulsive part---r at minimum well in the potential plot


saltConc = 0.6                   # salt concentration 0.6 moles/ltr
kappa = round(0.33*np.sqrt(saltConc)*sigma,3)    # debye screening length in reduced unit---for dielectric 80 and at temp 25 degree Celsius, kappa = 0.33(sqrt(I)) in per Angstrom unit
rcoul = 10.0/kappa               # to calculate cutoff for coulomb interaction, usually, kappa*cutoff = 10.






'''for incorporating hydrogen bonding, I tried to modify gauss potential
of the form     E = -A exp(-Br^2)
A is in unit of energy, B in unit of per distance square and we need cutoff... I choose it to be 2 sigma
For hydrogen bonding, energy would be 2-12 kBT or 6-30 KJpermole
I choose A to be 5.68kBT and B to be 0.1/sigma square '''




gaussA = 3.408/epsilon            # 1 KBT is 0.6kcal/mole---For now I choose 5kBT at the point where there is cutoff for excluded volume
gaussB = 0.1                    # 0.1/sigma square in reduced unit will be 0.1
gaussCut = (15*sigma)/sigma




'''time conversion'''
timeConversion = 4200*1000/(10**(-10))**2                   # conversion from kcal/Angstrom^2 to gram meter^2 /sec^2 meter^2
ljUnitTime  =  ((epsilon/(massReal* (sigma)**2)*timeConversion)**(1.0/2.0))/10**15  # conversion factor for time is 0.0003742  (epsilon/(mass* (sigma)^2))^(1/2) in unit per second 
tau = 10000 * ljUnitTime                                          # tau is 10 picosecod --the value mentioned here is in femtosecond unit  1 pico second = 1000 femtosecond
dt = round(0.0001*tau,5)


'''charge conversion'''
pivalue = np.pi  
permittivity = 8.85*10**(-12)      # Its unit is (C^2 sec^2)/(kg meter^3)
avogadoNum = 6.023 * 10**(23)
conversion = 4200*10**(-10)
ljUnitcharge = np.sqrt(4*pivalue*epsilon*sigma*permittivity*conversion/avogadoNum)     #charge conversion (4 pi permittivity sigma epsilon)^(1/2)
chargeValue = round((1.6*10**(-19))/ljUnitcharge, 2)                                            # value of 1 e charge in reduced unit 



kBT = 0.6                   # Value in kcal/mole in 300 K temperature 
tempReduced = kBT/epsilon
tempRelax = 0.2/epsilon







'''Define parameters needed for the data file'''

numResidue = 72
nResidue = 3*numResidue

#nTail = 30 

nBonds = 3*numResidue -1  #+nTail-1
nAngle = 3*numResidue -2 # +nTail-2          # Each residue have 2 backbone angles and one side angle ---first one has 2 and last one has 2
#nDihedral= 2*(numResidue-1) -1   # Each residue have 2 backbone dihedral(Phos-Sugar-Phos-Sugar, Sugar-Phos-Sugar-Phos) last residue doesnot start with, second last has only PSPS     
nDihedral= 0    
nImproper = 0


nAtomType = 25+1
nBondType = 1
nangleType = 2
#ndiheType = 1
ndiheType = 0
improType = 0




'''Structure'''

'''Read sequence file'''

inseq = open('sequence','r')
firstline = inseq.readline()
inseq.readline()

'''save the information as resNum, resType, charge, x, y, z'''

x = -0.5
y = 0.0
z = 68.

resInfo=[]
counter = 0

for index, value in enumerate(inseq):
    everyline = value.split()
    
    
    resNum = index+1
    #Phosphate
    counter = counter+1
    resType = "1"
    charge = -1.*chargeValue
    phos = [resNum, resType, charge, x, y, z-(index*1.13)]
    resInfo.append(phos)
    
    #Sugar
    counter = counter+1
    resType = "2"
    charge = 0.0*chargeValue
    sugar = [resNum, resType, charge, x+0.7, y, z+0.5-(index*1.13)]
    resInfo.append(sugar)
    
    #Bases
    counter = counter+1
    base= [resNum, int(everyline[0]), charge, x+1.5, y, z+0.5-(index*1.13)]
    resInfo.append(base)
    

nNucleic= counter




#'''get coordinates for tail atoms'''
#tailInfo = []
#resNum = numResidue+1
#x = 0.0
#y = 0.0
#z = 15.
#
#for index in range(nTail):
#    resType = "1"
#    charge = -1.0*chargeValue
#    tailCoord = [resNum, resType, charge, x, y, z-index-1]
#    tailInfo.append(tailCoord)
#








'''Pore coordinates---
 dimensions'''
 
resNum=resNum+1
charge =0.0*chargeValue

radius = 10.0/(2.0*sigma)   # Diameter of beta barrel in ahl is 1.4nm
height = 50.0/sigma         # Height of pore is choosen to be 5nm 
zvalue = int(height/2.)+1


'''I choose x and y range to go from -20 sigma unit to +20 sigma unit'''
xStart = -20
xEnd = 20
yStart = -20
yEnd = 20

'''(0,0,0) lies in the Up wall'''

zUp = height/2. 
zDown = -height/2.
dz = (zUp-zDown)*sigma             #in real unit

resType = 26
poreResType = resType


'''UpWall coordinates'''

UpWall = []
countUp = 0
for i in range(xStart,xEnd+1):
    for j in range(yStart,yEnd+1):
        value = [resNum, resType, charge, i, j, zUp]
        if (i**2 + j**2)>=(radius+1)**2:            # since center is (0, 0, z) and we don't care for z values for now (if (x-xc)^2 + (y -yc)^2 < radius) don't write the coordinates
            UpWall.append(value)
            countUp = countUp +1




'''Downwall coordinates'''
DownWall = []
countDown = 0
for i in range(xStart,xEnd+1):
    for j in range(yStart,yEnd+1):
        value = [resNum, resType, charge, i, j, zDown]
        if (i**2 + j**2)>=(radius+1)**2:
            DownWall.append(value)
            countDown = countDown +1


'''cylinder coordinates'''
countCylinder = 0
pore = []
for j in range(-int(height/2.),zvalue):
    for i in range(0,360,20):
        x = radius*np.cos(i*np.pi/180.)
        y = radius*np.sin(i*np.pi/180.)
        z = j
        value = [resNum, resType, charge, x, y, z]
        pore.append(value)
        countCylinder = countCylinder +1
    

nPore = countUp + countDown + countCylinder







#nTotal = nResidue+nTail+nPore
#nMobile = nResidue+nTail


nTotal = nResidue+nPore
nMobile = nResidue





'''Generate data file for the translocation simulation'''

outputatoms = open('system.data','w')
outNoPore = open('mobile.data','w')


outputatoms.write("LAMMPS Description\n\n")

outputatoms.write('{0} {1}\n'.format(nTotal, "atoms"))
outputatoms.write('{0} {1}\n'.format(int(nBonds), "bonds"))
outputatoms.write('{0} {1}\n'.format(int(nAngle), "angles"))
outputatoms.write('{0} {1}\n'.format(int(nDihedral), "dihedrals"))
outputatoms.write('{0} {1}\n\n'.format(int(nImproper), "impropers"))

outputatoms.write('{0} {1} {2}\n'.format(int(nAtomType), "atom", "types"))
outputatoms.write('{0} {1} {2}\n'.format(int(nBondType), "bond", "types"))
outputatoms.write('{0} {1} {2}\n'.format(int(nangleType), "angle", "types"))
outputatoms.write('{0} {1} {2}\n'.format(int(ndiheType), "dihedral", "types"))
outputatoms.write('{0} {1} {2}\n\n'.format(int(improType), "improper", "types"))


outputatoms.write('{0} {1} {2} {3}\n'.format("-36.0", "36.0", "xlo", "xhi" ))
outputatoms.write('{0} {1} {2} {3}\n'.format("-36.0", "36.0", "ylo", "yhi" ))
outputatoms.write('{0} {1} {2} {3}\n\n'.format("-140.0", "140.0", "zlo", "zhi" ))



outputatoms.write("Masses\n\n")
for j in range(nAtomType):
    outputatoms.write('{0} {1}\n'.format(int(j+1), massReduced))
    







outNoPore.write("LAMMPS Description\n\n")

outNoPore.write('{0} {1}\n'.format(nMobile, "atoms"))
outNoPore.write('{0} {1}\n'.format(int(nBonds), "bonds"))
outNoPore.write('{0} {1}\n'.format(int(nAngle), "angles"))
outNoPore.write('{0} {1}\n'.format(int(nDihedral), "dihedrals"))
outNoPore.write('{0} {1}\n\n'.format(int(nImproper), "impropers"))

outNoPore.write('{0} {1} {2}\n'.format(int(nAtomType), "atom", "types"))
outNoPore.write('{0} {1} {2}\n'.format(int(nBondType), "bond", "types"))
outNoPore.write('{0} {1} {2}\n'.format(int(nangleType), "angle", "types"))
outNoPore.write('{0} {1} {2}\n'.format(int(ndiheType), "dihedral", "types"))
outNoPore.write('{0} {1} {2}\n\n'.format(int(improType), "improper", "types"))


outNoPore.write('{0} {1} {2} {3}\n'.format("-36.0", "36.0", "xlo", "xhi" ))
outNoPore.write('{0} {1} {2} {3}\n'.format("-36.0", "36.0", "ylo", "yhi" ))
outNoPore.write('{0} {1} {2} {3}\n\n'.format("-140.0", "140.0", "zlo", "zhi" ))



outNoPore.write("Masses\n\n")
for j in range(nAtomType):
    outNoPore.write('{0} {1}\n'.format(int(j+1), massReduced))
    
    



   
outputatoms.write("\nAtoms\n\n")
outNoPore.write("\nAtoms\n\n")




'''Generate atom coordinates'''

'''coordinates of nucleotides'''



counter = 0
for index in range(nNucleic):
    counter = counter+1 
    outputatoms.write('{} {} {} {} {} {} {}\n'.format( int(counter), int(resInfo[index][0]), int(resInfo[index][1]), float(resInfo[index][2]), float(resInfo[index][3]), float(resInfo[index][4]), float(resInfo[index][5]) ))
    outNoPore.write('{} {} {} {} {} {} {}\n'.format( int(counter), int(resInfo[index][0]), int(resInfo[index][1]), float(resInfo[index][2]), float(resInfo[index][3]), float(resInfo[index][4]), float(resInfo[index][5]) ))




#for index in range(nTail):
#    counter = counter+1
#    outputatoms.write("{} {} {} {} {} {} {}\n".format( int(counter), int(tailInfo[index][0]), int(tailInfo[index][1]), float(tailInfo[index][2]), float(tailInfo[index][3]), float(tailInfo[index][4]), float(tailInfo[index][5]), ))
#    outNoPore.write("{} {} {} {} {} {} {}\n".format( int(counter), int(tailInfo[index][0]), int(tailInfo[index][1]), float(tailInfo[index][2]), float(tailInfo[index][3]), float(tailInfo[index][4]), float(tailInfo[index][5]), ))
    

for index in range(countUp):
    counter=counter+1
    outputatoms.write('{} {} {} {:0.3f} {} {} {}\n'.format( int(counter), int(UpWall[index][0]),  int(UpWall[index][1]), float(UpWall[index][2]), float(UpWall[index][3]), float(UpWall[index][4]), UpWall[index][5]))
   
   
for index in range(countDown):
    counter=counter+1
    outputatoms.write('{} {} {} {:0.3f} {} {} {}\n'.format( int(counter), int(DownWall[index][0]),  int(DownWall[index][1]), float(DownWall[index][2]), float(DownWall[index][3]), float(DownWall[index][4]), DownWall[index][5]))
   
   
for index in range(countCylinder):
    counter=counter+1
    outputatoms.write('{} {} {} {:0.3f} {} {} {}\n'.format( int(counter), int(pore[index][0]), int(pore[index][1]), float(pore[index][2]),  float(pore[index][3]), float(pore[index][4]), float(pore[index][5])))







poreAtomTypelett = "C"
outputpore = open("pore.xyz", "w")
outputpore.write("{}\n".format(nPore))
outputpore.write("Atoms\n")
for index in range(countUp):
    outputpore.write('{} {} {} {}\n'.format(poreAtomTypelett, float(UpWall[index][3]), float(UpWall[index][4]), UpWall[index][5]))
for index in range(countDown):
    outputpore.write('{} {} {} {}\n'.format(poreAtomTypelett, float(DownWall[index][3]), float(DownWall[index][4]), DownWall[index][5]))
for index in range(countCylinder):
    outputpore.write('{} {} {} {}\n'.format(poreAtomTypelett, float(pore[index][3]), float(pore[index][4]), float(pore[index][5])))










    

'''Generate list of bonds'''

outputatoms.write('\n Bonds \n\n')
outNoPore.write('\n Bonds \n\n')

bondType = "1"  
num = 0
for index in range(numResidue-1):
    num = num +1
    outputatoms.write('{0} {1} {2} {3}\n'.format(int(num), bondType, int((index+1)*3 - 2), int(index+1)*3 -1))
    outNoPore.write('{0} {1} {2} {3}\n'.format(int(num), bondType, int((index+1)*3 - 2), int(index+1)*3 -1))
    num = num +1
    outputatoms.write('{0} {1} {2} {3}\n'.format(int(num), bondType, int((index+1)*3 - 1), int(index+1)*3 -0))
    outNoPore.write('{0} {1} {2} {3}\n'.format(int(num), bondType, int((index+1)*3 - 1), int(index+1)*3 -0))
    num = num+1
    outputatoms.write('{0} {1} {2} {3}\n'.format(int(num), bondType, int((index+1)*3 - 1), int(index+1)*3 +1))
    outNoPore.write('{0} {1} {2} {3}\n'.format(int(num), bondType, int((index+1)*3 - 1), int(index+1)*3 +1))


'''last bond'''
index =  int((index+1)*3 +1)
num = num+1
outputatoms.write('{0} {1} {2} {3}\n'.format(int(num), bondType, int(index), int(index+1) ))
outNoPore.write('{0} {1} {2} {3}\n'.format(int(num), bondType, int(index), int(index+1) ))
num = num +1
outputatoms.write('{0} {1} {2} {3}\n'.format(int(num), bondType, int(index+1), int(index+2) ))
outNoPore.write('{0} {1} {2} {3}\n'.format(int(num), bondType,int(index+1), int(index+2) ))




#'''bonds attached to the tail'''
#for i in range(num,num+nTail-1):
#    outputatoms.write('{0} {1} {2} {3}\n'.format(int(i+1), bondType, int(i+1), int(i +2)))
#    outNoPore.write('{0} {1} {2} {3}\n'.format(int(i+1), bondType, int(i+1), int(i +2)))







'''Generate list of angles'''

outputatoms.write("\nAngles\n\n")
outNoPore.write("\nAngles\n\n")
num =0

'''backbone angles  --Phosphate-Sugar-Phospate, Sugar-Phosphate-Sugar'''

angleType = "1"

'''first angle'''
num =num+1
outputatoms.write('{} {} {} {} {}\n'.format( int(num), angleType, "1", "2", "4"))
outNoPore.write('{} {} {} {} {}\n'.format( int(num), angleType, "1", "2", "4"))
'''in between'''
for index in range(numResidue-2):
    num =num+1
    outputatoms.write('{} {} {} {} {}\n'.format( int(num), angleType, int((index+1)*3 - 1), int((index+1)*3 +1), int((index+1)*3+2) ))
    outNoPore.write('{} {} {} {} {}\n'.format( int(num), angleType, int((index+1)*3 - 1), int((index+1)*3 +1), int((index+1)*3+2) ))
    num = num+1
    outputatoms.write('{} {} {} {} {}\n'.format( int(num), angleType, int((index+1)*3 + 1), int((index+1)*3 +2), int((index+1)*3+4) ))
    outNoPore.write('{} {} {} {} {}\n'.format( int(num), angleType, int((index+1)*3 + 1), int((index+1)*3 +2), int((index+1)*3+4) ))



'''last angle'''

index = int((index+1)*3 +4)
num = num+1
outputatoms.write('{} {} {} {} {}\n'.format( int(num), angleType, int(index-2), int(index), int(index+1) ))
outNoPore.write('{} {} {} {} {}\n'.format( int(num), angleType, int(index-2), int(index), int(index+1) ))
#
#'''angle of the attached tail'''
#count= index+1
#
#for i in range(nTail-2):
#    num=num+1
#    outputatoms.write('{} {} {} {} {}\n'.format( int(num), angleType, int(count+i-1), int(count+i), int(count+i+1)  ))
#    outNoPore.write('{} {} {} {} {}\n'.format( int(num), angleType, int(count+i-1), int(count+i), int(count+i+1)  ))




'''side angles  ---Phos-Sugar-Base for each nucleotide'''

angleType ="2"
for index in range(numResidue):
    num = num +1
    outputatoms.write('{0} {1} {2} {3} {4}\n'.format(int(num), angleType, int((index+1)*3 - 2), int((index+1)*3 -1), int((index+1)*3-0)))
    outNoPore.write('{0} {1} {2} {3} {4}\n'.format(int(num), angleType, int((index+1)*3 - 2), int((index+1)*3 -1), int((index+1)*3-0)))




#'''Generate dihedral angles'''
#'''Backbone ---PSPS and SPSP'''
#
#outputatoms.write("\n Dihedrals\n\n")
#outNoPore.write("\n Dihedrals\n\n")
#num =0
#diheType="1"
#for index in range(numResidue -2):
#    num = num+1
#    outputatoms.write('{} {} {} {} {} {}\n'.format( int(num), diheType, int((index+1)*3 -2), int((index+1)*3 -1), int((index+1)*3 +1), int((index+1)*3 +2) ))
#    outNoPore.write('{} {} {} {} {} {}\n'.format( int(num), diheType, int((index+1)*3 -2), int((index+1)*3 -1), int((index+1)*3 +1), int((index+1)*3 +2) ))
#    num = num+1
#    outputatoms.write('{} {} {} {} {} {}\n'.format( int(num), diheType, int((index+1)*3 -1), int((index+1)*3 +1), int((index+1)*3 +2), int((index+1)*3 +4) ))
#    outNoPore.write('{} {} {} {} {} {}\n'.format( int(num), diheType, int((index+1)*3 -1), int((index+1)*3 +1), int((index+1)*3 +2), int((index+1)*3 +4) ))
#
#'''last one'''
#
#num = num+1
#index=index+1
#outputatoms.write('{} {} {} {} {} {}\n'.format( int(num), diheType, int((index+1)*3 -2), int((index+1)*3 -1), int((index+1)*3 +1), int((index+1)*3 +2) ))
#outNoPore.write('{} {} {} {} {} {}\n'.format( int(num), diheType, int((index+1)*3 -2), int((index+1)*3 -1), int((index+1)*3 +1), int((index+1)*3 +2) ))





























'''Generate inputfile for the simulation'''
repeatline = "##---------------------------------\n"
outInFile = open('system.in','w')
outInFile.write("################ RNA with tail ##########\n")
outInFile.write("clear\n\n")

outInFile.write("##Initialization\n")
outInFile.write(repeatline)

outInFile.write('{0} {1}\n'.format("units", "lj"))
outInFile.write('{0} {1}\n'.format("dimension", "3"))
outInFile.write('{0} {1}\n'.format("atom_style", "full"))
outInFile.write('{0} {1} {2} {3}\n\n'.format("boundary", "s", "s", "s"))

outInFile.write("##Atom definition\n")
outInFile.write(repeatline)
outInFile.write('{0} {1}\n\n'.format("read_data", "system.data"))

outInFile.write("## Force definition\n")
outInFile.write(repeatline)

outInFile.write("## Bond definition\n")
outInFile.write('{0} {1}\n'.format("bond_style", "harmonic"))
outInFile.write('{0} {1} {2} {3}\n\n'.format("bond_coeff", bondType, kb, l0))

outInFile.write("## Angle definition\n")
outInFile.write('{0} {1}\n'.format("angle_style", "harmonic"))
outInFile.write('{0} {1} {2} {3}\n'.format("angle_coeff", bondType, kTheta, theta0))
outInFile.write('{0} {1} {2} {3}\n\n'.format("angle_coeff", angleType, kTheta, theta0))

#outInFile.write("## Dihedral definition\n")
#outInFile.write('{0} {1}\n'.format("dihedral_style", "helix"))
#outInFile.write('{0} {1} {2} {3} {4}\n\n'.format("dihedral_coeff", diheType, diheA, diheB, diheC))


outInFile.write("special_bonds	lj/coul 0 1 1 angle yes dihedral yes\n")
outInFile.write('{0} {1}\n'.format("dielectric", dielectric))

outInFile.write("## Pair definition\n")
outInFile.write(repeatline)

outInFile.write('{0} {1} {2} {3} {4} {5} {6:0.3f} {7} {8}\n'.format("pair_style", "hybrid/overlay", "lj/cut", rlj, "coul/debye", kappa, rcoul, "gauss", gaussCut ))

outInFile.write("pair_modify	shift yes\n")

'''pair coefficient for lj'''

for i in range(nAtomType):
    for j in range(nAtomType):
        if i<j:
            break
        else:
            
            if i==j and i>=3:
                outInFile.write('{0} {1} {2} {3} {4} {5} {6:0.3f}\n'.format("pair_coeff", j+1, i+1, "lj/cut", float(epsilonlj/epsilon), float(sigmalj/sigma), rlocal ))
                outInFile.write("pair_coeff {} {} gauss {} {}\n".format( j+1, i+1, gaussA, gaussB))
            else:
                outInFile.write('{0} {1} {2} {3} {4} {5} {6:0.3f}\n'.format("pair_coeff", j+1, i+1, "lj/cut", float(epsilonlj/epsilon), float(sigmalj/sigma), rlocal ))

'''pair coefficient for coulomb/debye'''
outInFile.write('{0} {1} {2} {3}\n'.format("pair_coeff", 1, 1, "coul/debye"))


run = 30000000
#run = 1000

endpartStr ="\n#Group Definition \n#---------------------------------------------------------\n"
endpartStr +="group		mobile molecule <> 1 {}\n".format(numResidue)
endpartStr +="group		pore molecule <> {} {}\n".format(numResidue+1, numResidue+1)
endpartStr += "group         initial id {}\n\n".format(nMobile-2)

endpartStr +="#Neighbor Modify \n#-------------------------------------------------------\n"
endpartStr +="neigh_modify	exclude type {} {}\n\n".format(poreResType, poreResType)


endpartStr += "\n#Timestep etc\n"
endpartStr += "#-----------------------------------------------\n"
endpartStr +="timestep                {}\n".format(dt)
endpartStr +="run_style               verlet\n"
endpartStr +="#velocity               mobile create {} 5748\n".format(tempRelax)
endpartStr += "velocity 	initial  zero  linear\n\n"


endpartStr+="#Fix\n#------------------------------------------\n"
endpartStr+="fix 2 mobile nve molecule\n"
endpartStr+="fix 3 all langevin {} {} 1.0 23541\n".format(tempRelax, tempRelax)
endpartStr+="fix	4 initial setforce 0.0 0.0 0.0\n\n"


endpartStr+="#Dump\n#----------------------------------------------\n"
endpartStr+="thermo_style            custom step temp evdwl ecoul ebond eangle pe ke etotal \n"
endpartStr+="thermo          1000\n"
endpartStr+="dump  1 mobile custom 5000 minimize.lammpstrj id mol type x y z \n"
endpartStr+="dump  5 all xyz {} dump.xyz\n".format(run)
#endpartStr+="dump  5 mobile xyz 1 dump.xyz\n"

endpartStr+="run     {}\n".format(run)
endpartStr+="#write_data last.data\n"

outInFile.write(endpartStr)





    
    

outputatoms.close()
outNoPore.close()
outputpore.close()
outInFile.close()
outputatoms.close()
