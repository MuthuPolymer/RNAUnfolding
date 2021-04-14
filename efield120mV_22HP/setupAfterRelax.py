# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 18:43:42 2018

@author: sana
"""

'''Generating input and data file to apply electric field to the system '''
'''lastData.xyz contains coordinates and atom type of the relaxed system'''



import numpy as np






'''conversion from real to lj unit'''

sigma = 2.5         #bead size in angstrom unit
epsilon = 0.2       # in kcal/mol

epsilonH = 1.0/epsilon    # strength of hydrophobic interaction        

massReal = 96.0
massReduced = 96.0/massReal                       # mass is 1.0 is reduced unit. It is taken each bead is of same mass i.e. 96 Da i.e. 96 gram/mole

#Bond coefficients

l0 = 3.0/sigma             # equilibrium bond length 3.0 angstrom unit
#kb = round((1./2.)*100.*epsilonH/l0**2.,3)     # 1/2 100 epsilonH /a^2
kb = 171.0 * (sigma**2)/epsilon   # harmonic constant for the harmonic potential. value is 171 kcal/mol/Angstrom^2 ( energy/distance^2)kb

#Angle coefficients

#kTheta = (1./2.)*20.*epsilonH      # in unit (energy/radian^2)
kTheta = 50./epsilon                # in unit of energy----Our angle potential is cosine/squared K(costheta - costheta0)^2
theta0 = 105.                       # in degrees



#Dihedral coefficients 
#
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




gaussA = 6.612/epsilon            # 1 KBT (20 degree C) is 0.58kcal/mole---For now I choose 10kBT at the point where there is cutoff for excluded volume
gaussB = 0.1                    # 0.1/sigma square in reduced unit will be 0.1
gaussCut = (15*sigma)/sigma



'''time conversion'''
timeConversion = 4200*1000/(10**(-10))**2                   # conversion from kcal/Angstrom^2 to gram meter^2 /sec^2 meter^2
ljUnitTime  =  ((epsilon/(massReal* (sigma)**2)*timeConversion)**(1.0/2.0))/10**15  # conversion factor for time is 0.0003742  (epsilon/(mass* (sigma)^2))^(1/2) in unit per second 
tau = 10000 * ljUnitTime                                          # tau is 10 picosecod --the value mentioned here is in femtosecond unit  1 pico second = 1000 femtosecond
dt = 0.00005*tau                                         # dt = 0.5 femtosecond




'''charge conversion'''
pivalue = np.pi  
permittivity = 8.85*10**(-12)      # Its unit is (C^2 sec^2)/(kg meter^3)
avogadoNum = 6.023 * 10**(23)
conversion = 4200*10**(-10)
ljUnitcharge = np.sqrt(4*pivalue*epsilon*sigma*permittivity*conversion/avogadoNum)     #charge conversion (4 pi permittivity sigma epsilon)^(1/2)
chargeValue = round((1.6*10**(-19))/ljUnitcharge, 3)                                            # value of 1 e charge in reduced unit 



kBT = 0.58                  # Value in kcal/mole in 293 K  i.e. 20 C temperature 
tempReduced = kBT/epsilon




'''Pore  dimensions'''


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



'''efield conversion'''
ljUnitEfield = np.sqrt(4*pivalue*permittivity*sigma*epsilon*conversion/avogadoNum)*(avogadoNum/4200.)*(sigma/epsilon)

efieldValue = ((120.0/1000.)/dz)*ljUnitEfield               #120 mV









'''Define parameters needed for the data file'''

numResidue = 72
nResidue = 3*numResidue

nMobile = nResidue



nBonds = 3*numResidue-1
nAngle = 3*numResidue -2          # Each residue have 2 backbone angles and one side angle ---first one has 2 and last one has 2
#nDihedral= 2*(numResidue-1) -1   # Each residue have 2 backbone dihedral(Phos-Sugar-Phos-Sugar, Sugar-Phos-Sugar-Phos) last residue doesnot start with, second last has only PSPS     
nDihedral= 0    
nImproper = 0


nAtomType = 25+1
nBondType = 1
nangleType = 2
#ndiheType = 1
ndiheType = 0
improType = 0
























'''Read data.xyz and get coordinates'''


readdata= open('data.xyz','r')
firstline = readdata.readline()
secondline = readdata.readline()
inputdata = readdata.readlines()   #coordinates data file

nTotal = int(firstline)  

'''Generate data file for the translocation simulation'''

outputatoms = open('efield.data','w')

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


outputatoms.write('{0} {1} {2} {3}\n'.format("-40.0", "40.0", "xlo", "xhi" ))
outputatoms.write('{0} {1} {2} {3}\n'.format("-40.0", "40.0", "ylo", "yhi" ))
outputatoms.write('{0} {1} {2} {3}\n\n'.format("-140.0", "140.0", "zlo", "zhi" ))


outputatoms.write("Masses\n\n")
for j in range(nAtomType):
    outputatoms.write('{0} {1}\n'.format(int(j+1), massReduced))

outputatoms.write("\nAtoms\n\n")

'''split contents of lines from third line of xyz file'''

chargePhos = -1.*chargeValue

coordAll=np.zeros((nTotal,6))
residueType=np.zeros(nTotal)
counter = 1

for i in range(nTotal):
    if i<int(nResidue):
        residueType[i]=int(round((i)/3) +1)
    else:
        residueType[i] = int(numResidue+1)
    
    

for idx,lines in enumerate(inputdata):
    everyline = lines.split()    
    if everyline[0]=="1":
        coordAll[idx][0]=counter
        coordAll[idx][1]=float(everyline[0])
        coordAll[idx][2]=round(chargePhos,3)
        coordAll[idx][3]=float(everyline[1])
        coordAll[idx][4]=float(everyline[2])
        coordAll[idx][5]=float(everyline[3])        
    else:
        coordAll[idx][0]=counter
        coordAll[idx][1]=float(everyline[0])
        coordAll[idx][2]=0.0
        coordAll[idx][3]=float(everyline[1])
        coordAll[idx][4]=float(everyline[2])
        coordAll[idx][5]=float(everyline[3])
    counter = counter+1

for index in range(nTotal):
    outputatoms.write("{} {} {} {} {} {} {}\n".format(int(coordAll[index][0]), int(residueType[index]), int(coordAll[index][1]), coordAll[index][2], coordAll[index][3], coordAll[index][4], coordAll[index][5])) 
   
    
'''Generate bonds and angles'''

outputatoms.write('\n Bonds \n\n')



bondType = "1"  
num = 0
for index in range(numResidue-1):
    num = num +1
    outputatoms.write('{0} {1} {2} {3}\n'.format(int(num), bondType, int((index+1)*3 - 2), int(index+1)*3 -1))
    num = num +1
    outputatoms.write('{0} {1} {2} {3}\n'.format(int(num), bondType, int((index+1)*3 - 1), int(index+1)*3 -0))
    num = num+1
    outputatoms.write('{0} {1} {2} {3}\n'.format(int(num), bondType, int((index+1)*3 - 1), int(index+1)*3 +1))

'''last bond'''
index =  int((index+1)*3 +1)
num = num+1
outputatoms.write('{0} {1} {2} {3}\n'.format(int(num), bondType, int(index), int(index+1) ))
num = num +1
outputatoms.write('{0} {1} {2} {3}\n'.format(int(num), bondType, int(index+1), int(index+2) ))






'''Generate list of angles'''

outputatoms.write("\nAngles\n\n")
num =0

'''backbone angles  --Phosphate-Sugar-Phospate, Sugar-Phosphate-Sugar'''

angleType = "1"

'''first angle'''
num =num+1
outputatoms.write('{} {} {} {} {}\n'.format( int(num), angleType, "1", "2", "4"))
'''in between'''
for index in range(numResidue-2):
    num =num+1
    outputatoms.write('{} {} {} {} {}\n'.format( int(num), angleType, int((index+1)*3 - 1), int((index+1)*3 +1), int((index+1)*3+2) ))
    num = num+1
    outputatoms.write('{} {} {} {} {}\n'.format( int(num), angleType, int((index+1)*3 + 1), int((index+1)*3 +2), int((index+1)*3+4) ))



'''last angle'''

index = int((index+1)*3 +4)
num = num+1
outputatoms.write('{} {} {} {} {}\n'.format( int(num), angleType, int(index-2), int(index), int(index+1) ))


'''side angles  ---Phos-Sugar-Base for each nucleotide'''

angleType ="2"
for index in range(numResidue):
    num = num +1
    outputatoms.write('{0} {1} {2} {3} {4}\n'.format(int(num), angleType, int((index+1)*3 - 2), int((index+1)*3 -1), int((index+1)*3-0)))




'''Generate in file'''


outInFile = open('efield.in','w')
'''Generate input file for the simulation'''

repeatline = "##---------------------------------\n"
outInFile.write("################ RNA with tail translocation under electric field ##########\n")
outInFile.write("clear\n\n")

outInFile.write("##Initialization\n")
outInFile.write(repeatline)

outInFile.write('{0} {1}\n'.format("units", "lj"))
outInFile.write('{0} {1}\n'.format("dimension", "3"))
outInFile.write('{0} {1}\n'.format("atom_style", "full"))
outInFile.write('{0} {1} {2} {3}\n\n'.format("boundary", "p", "p", "f"))


outInFile.write("##Atom definition\n")
outInFile.write(repeatline)
outInFile.write('{0} {1}\n\n'.format("read_data", "efield.data"))


outInFile.write("## Force definition\n")
outInFile.write(repeatline)

outInFile.write("## Bond definition\n")
outInFile.write('{0} {1}\n'.format("bond_style", "harmonic"))
outInFile.write('{0} {1} {2} {3}\n\n'.format("bond_coeff", bondType, kb, l0))

outInFile.write("## Angle definition\n")
outInFile.write('{0} {1}\n'.format("angle_style", "cosine/squared"))
outInFile.write('{0} {1} {2} {3}\n'.format("angle_coeff", bondType, kTheta, theta0))
outInFile.write('{0} {1} {2} {3}\n\n'.format("angle_coeff", angleType, kTheta, theta0))

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


run = 50000000

endpartStr ="\n#Group Definition \n#---------------------------------------------------------\n"
endpartStr +="group		mobile molecule <> 1 {}\n".format(numResidue)
endpartStr +="group		pore molecule <> {} {}\n".format(numResidue+1, numResidue+1)


endpartStr +="#Neighbor Modify \n#-------------------------------------------------------\n"
endpartStr +="neigh_modify	exclude type {} {}\n\n".format(poreResType, poreResType)


endpartStr += "\n#Timestep etc\n"
endpartStr += "#-----------------------------------------------\n"
endpartStr +="timestep                {}\n".format(dt)
endpartStr +="run_style               verlet\n"
endpartStr +="velocity                mobile create {} 65748\n\n".format(tempReduced)

endpartStr +="#Variable Definitions \n#----------------------------------------------\n"
endpartStr+="variable  y1 equal z[1]\n"
endpartStr+="variable  y2 equal z[4]\n"
endpartStr+="variable  y3 equal z[7]\n"
endpartStr+="variable  y4 equal z[10]\n"
endpartStr+="variable  y5 equal z[13]\n"
endpartStr+="variable  y6 equal z[16]\n"
endpartStr+="variable  y7 equal z[19]\n"
endpartStr+="variable  y8 equal z[22]\n"
endpartStr+="variable  y9 equal z[25]\n"
endpartStr+="variable  y10 equal z[28]\n"
endpartStr+="variable  y11 equal z[{}]\n\n".format(nMobile)

endpartStr += "variable  b1 equal {}\n".format(zDown)
endpartStr += "variable  b2 equal {}\n\n".format(zUp+5.)

endpartStr += "variable  e1 equal 0.0\n"
endpartStr += "variable  e2 equal 0.0\n"
endpartStr += "variable  e3 equal {}\n\n".format(efieldValue)

#endpartStr +="variable field atom ((v_e1*(z<{}))+(v_e2*(z>{}))+(v_e3*((x<{})&&(x>{})&&(y<{})&&(y>{})&&(z<{})&&(z>{}))))\n\n".format(zDown, zUp, -radius, radius, -radius, radius, zDown, zUp )
endpartStr +="variable field atom ((v_e3*((x>{})&&(x<{})&&(y>{})&&(y<{})&&(z>{})&&(z<{}))))\n\n".format( -radius, radius, -radius, radius, zDown, zUp )


endpartStr+="#Fix\n#------------------------------------------\n"
endpartStr+="fix 2 mobile nve molecule\n"
endpartStr+="fix 3 all langevin {} {} 1.0 23541\n".format(tempReduced, tempReduced)
endpartStr+="fix ef mobile efield 0.0 0.0 v_field\n\n"


endpartStr+="#Computes\n#----------------------------------------------\n"
endpartStr+="compute  rr  mobile  gyration\n"

endpartStr+="#Dump\n#----------------------------------------------\n"
endpartStr+="thermo_style            custom step temp evdwl ecoul ebond eangle pe ke etotal c_rr \n"
endpartStr+="thermo          1000\n"
endpartStr+="dump  1 mobile custom 5000 translocate.lammpstrj id mol type x y z \n"
endpartStr+="dump  5 mobile xyz 1000 dump.xyz\n\n"

endpartStr+="run  {}  ".format(run)
endpartStr+="every 1000 \"if '${y1}<${b1} && ${y2}<${b1} && ${y3}<${b1} && ${y4}<${b1} && ${y5}<${b1} && ${y6}<${b1} && ${y7}<${b1} && ${y8}<${b1} && ${y9}<${b1} && ${y10}<${b1} || ${y11}>${b2} ' then quit\"\n\n"
endpartStr+="#----------End-----------\n"

outInFile.write(endpartStr)










readdata.close()
outputatoms.close()
outInFile.close()
