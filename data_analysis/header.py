
#############################################3
#
#notes: code wil panic is partcile number drops to 0 beocuse of the len(i) in line 51
#
##############################################


import numpy as np
from scipy import constants

import matplotlib.pyplot as plt
from timeit import timeit


def get_coordinates(dat):
    a = dat.index(',')
    location = [dat.index('('), a, dat.index(',', a+1), dat.index(')')]
    coordinates = [ float(dat[(location[i]+1):location[i+1]]) for i in range(3) ]
    return coordinates


def get_particle_data( file ):

    Step, ParticleNumber, ParticleData = [], [], []

    with open( file , 'r') as myfile:

        while( True ):

            #Getting file data line-by-line.
            MetaData = myfile.readline()
            #Check end of file.
            if len(MetaData) == 0:
                break

            Step_tmp           = int(MetaData[ MetaData.index('-')+1 : MetaData.index(',') ])
            ParticleNumber_tmp = int(MetaData[ MetaData.index(',')+1 :  ])

            ParticleID       = []
            ParticleLocation_tmp = []

            for i in range(ParticleNumber_tmp):
                PositionData = myfile.readline()
                ParticleID.append(PositionData[ PositionData.index(',')+1: PositionData.index(':') ])
                ParticleLocation_tmp.append(get_coordinates(PositionData[ PositionData.index(' ')+1: PositionData.index('\n') ]))

            Step.append(Step_tmp)
            ParticleNumber.append(ParticleNumber_tmp)
            ParticleData.append(ParticleLocation_tmp)

    return Step, ParticleNumber, ParticleData

def get_Vrms( ParticleVelovity ):
    '''Gets root mean squeard velcity per run
    '''
    Vrms = [ np.sqrt(sum( [j[0]**2 + j[1]**2 + j[2]**2 for j in i] )/len(i)) if len(i) > 0 else 0.0 for i in ParticleVelovity]
    return Vrms

def get_kinetic_energy( ParticleVelovity, Mass):
    '''Gets kinetic_energy using root mean squeard velcity per run
    '''
    Vrms = get_Vrms(ParticleVelovity)
    KineticEnergy = [0.5*Mass*(v)**2 for v in Vrms]
    return KineticEnergy

def get_tempeture(ParticleVelovity, Mass):
    '''Gets temp[K] using kinetic_energy using root mean squeard velcity per run
    '''
    kB = constants.value(u'Boltzmann constant') #1.380649e-23 J K^-1
    KineticEnergy = get_kinetic_energy( ParticleVelovity, Mass)
    T = [ (2/3)*(Ke/kB) for Ke in KineticEnergy ]
    return T

if __name__ == "__main__":
    # s = """get_particle_data("../vel.txt")"""
    # print(timeit(stmt=s, number=10000, setup="from __main__ import get_particle_data"))
    fileV = "../vel_dipole.txt"
    fileP = "../pos.txt"

    VelData = get_particle_data( fileV )[2]
    print(get_tempeture( VelData, 87*1.66e-27 ))

    plt.plot(get_particle_data( fileV )[0], get_tempeture( VelData, 87*1.66e-27 ))
    plt.show()
