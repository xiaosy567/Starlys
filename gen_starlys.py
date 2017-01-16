import math
import numpy as np
from math import cos, sin, sqrt
from numpy import matrix
from random import randint
from numpy import array
from numpy import cross
from numpy import dot
from numpy import abs
import re
import string

import sys
sys.path.append('./') 
sys.path.append('./mymath') 
import PDB_module
import GRO_TOP
import ROT_AAL

def pointsOnSphere(N):
    N = float(N) # in case we got an int which we surely got
    pts = []
 
    inc = math.pi * (3 - math.sqrt(5))
    off = 2 / N
    for k in range(0, int(N)):
        y = k * off - 1 + (off / 2)
        r = math.sqrt(1 - y*y)
        phi = k * inc
        pts.append([math.cos(phi)*r, y, math.sin(phi)*r])
 
    return pts

def getX(x,scale,shift=0):
    return x*scale+shift

class DnaNanoParticle:
    '''
    This class is going to define a dna grafted particle
    '''

    def __init__(self, DnpAtomList=[], DnpBondList=[], DnpAngleList=[], 
        NpBall = 40, Ndna=2, RadiurBall = 8.0, BondSP = 1.0, pdbname = "1.pdb", fatom1 = 1, fatom2 = 20, itpname = "1.itp"):
	if(len(DnpAtomList) >=1 or len(DnpBondList)>=1 or len(DnpAngleList)>=1):
            print "----------------------------------------------------"
            print "The initial DNA nano particle class should be ENPTY!"
            print "----------------------------------------------------"
        self.DnpAtomList=[]
        self.DnpAtomITPList=[]
        self.DnpBondList=[]
        self.DnpAngleList=[]
        self.DnpConstraintsList=[]
        self.DnpdihedralList=[]
        #number of atom of a big ball
        N=NpBall
        #number of mao per big ball
        L=Ndna
        #radius of big ball
        s=RadiurBall
        ssp=BondSP
        lensstrand=1
        #reference atom index:
        refindex1=fatom1
        refindex2=fatom2
      
        mol = []
        bnd = []
        ang = []
        ndx = 0

        ##topology items
        moleculetype=[]
        atoms=[]
        bonds=[]
        angles=[]
        constraints=[]
        dihedrals=[]

        #NUMAtomType---> numform atom type
        NUMRIDUEType=0
     


        ##read pdb
        atomlist=GRO_TOP.Read_pdb(pdbname)
        print "+++ pdb info +++"
        print atomlist
        print "+++ pdb info +++"

        ##read_itp
        toplist=GRO_TOP.Read_itp(itpname)
        print "+++ itp info +++"
        print toplist
        print "+++ itp info +++"



        #####################
        #generate a big ball#
        #####################
        #RESIDUE TYPE
        NUMRIDUEType = 1

        for n,pt in enumerate(pointsOnSphere(N)):
            ndx = ndx+1
            #mol formate (ndx, numform-atom-type, string-form-of-atom-type, residue-type, x, y z)
            mol.append([ndx, NUMRIDUEType, 'NP', 'DNP', getX(pt[0],s), getX(pt[1],s), getX(pt[2],s)])
            #add atom stup for itp
            #; nr type  resnr residu  atom  cgnr charge
            atoms.append([ndx, 'NP', 1, 'BAL', 'NN', ndx,  0.0]) 
        
        ###################################
        #generate the grafted dna branches#
        ###################################
        #s = s+sss
        roots = pointsOnSphere(L)
        # nroots = []
        for n,pt in enumerate(roots):
            
            print "..................add ", n, "th ...................."
            

            ####################
            #generate positon for mao attachment
            
            pt_pos=np.array([getX(pt[0],s+ssp), getX(pt[1],s+ssp), getX(pt[2],s+ssp)])

            #GET the neareast particle on the surface of big ball
            p1=array([pt[0], pt[1], pt[2]])
            pointer=0
            pCos=-1
            for m in range(N):
                p2=array([mol[m][4], mol[m][5], mol[m][6]])
                if(dot(p1, p2)> pCos):
                    pCos=dot(p1, p2)
                    pointer=m

            ########################################################
            #add mao
            ########################################################
            #axis (origin)->(particle on surface)
            axis_a1 = np.array([0, 0, 0])
            axis_a2 = pt_pos
            axis_a = axis_a2 - axis_a1
            # print "axis_a1:", axis_a1
            # print "axis_a2:", axis_a2
            # print "axis_a:", axis_a

            #ref of mao (a1)->(a2)
            ref_a1 = np.array([atomlist[refindex1][5], atomlist[refindex1][6], atomlist[refindex1][7] ])
            ref_a2 = np.array([atomlist[refindex2][5], atomlist[refindex2][6], atomlist[refindex2][7] ])
            ref_a = ref_a2 - ref_a1
            # print "ref_a1:", ref_a1
            # print "ref_a2:", ref_a2
            # print "ref_a:", ref_a
            
            #rotation angle
            rotate_angle = ROT_AAL.angle_vectors(axis_a, ref_a)

            #rotation axis (ref_a X axis_a)
            axis_v = ROT_AAL.return_unitvector(np.cross(ref_a, axis_a))
            
            #the distance between ref_a1 and pt_pos
            transmove = ref_a1 - pt_pos
            # print "transmove:", transmove
            # print "axis_v:", axis_v
            # print "rotate_angle:", rotate_angle

            # zerovector = np.array([0, 0, 0])
            # ref_a1_array = np.array([atomlist[refindex1][5], atomlist[refindex1][6], atomlist[refindex1][7], 1.0 ])
            # ref_a1_m = ROT_AAL.rotate(zerovector, axis_v, ref_a1_array, rotate_angle)
            # transmove = np.array([ref_a1_m[0]-pt_pos[0], ref_a1_m[1]-pt_pos[1], ref_a1_m[2]-pt_pos[2] ])
            
            ndx_atomlist = ndx # atom index in the atomlist
            typeindex_atomlist = NUMRIDUEType
            for j, pm in enumerate(atomlist):
                atomindex = ndx_atomlist + 1 + j
                RESIUDETYPE_INDEX = typeindex_atomlist + pm[4] #residue index in the molecule list
                
                # move the system to the point of axis_a
                pmarray=np.array([ pm[5] - transmove[0], pm[6] - transmove[1], pm[7] - transmove[2], 1.0 ]);
                pam = ROT_AAL.rotate(pt_pos, axis_v + pt_pos, pmarray, rotate_angle)
                mol.append([atomindex, RESIUDETYPE_INDEX, pm[2], pm[3], pam[0], pam[1], pam[2] ])
                if j==refindex1:
                    bonds.append([pointer+1, atomindex, 1, 0.33, 5000.0]);
                elif j == (len(atomlist)-1):
                    NUMRIDUEType = NUMRIDUEType + pm[4]
                    ndx = ndx + len(atomlist)

            ##generate itp
            for j, p in enumerate(toplist):
                name=p[0];
                elements=p[1];
                print '[ ', name, ' ]';
                if name == 'atoms':
                    for ele in elements:
                        #1    Qd     1   LYS    BB     1  1.0000
                        atoms.append([ ele[0] + ndx_atomlist, ele[1], ele[2]+typeindex_atomlist , ele[3], ele[4], ele[5] + ndx_atomlist, ele[6] ]);
                        print '\t', ele[0] + ndx_atomlist, '\t', ele[1], '\t', ele[2], '\t', ele[3], '\t', ele[4], '\t', ele[5], '\t', ele[6];
                elif name == 'bonds':
                    for ele in elements:
                        bonds.append([ ele[0] + ndx_atomlist, ele[1] + ndx_atomlist, ele[2], ele[3], ele[4] ]);
                        print '\t', ele[0] + ndx_atomlist, '\t', ele[1] + ndx_atomlist, '\t', ele[2], '\t', ele[3], '\t', ele[4];
                elif name == 'constraints':
                    for ele in elements:
                        constraints.append([ ele[0] + ndx_atomlist, ele[1] + ndx_atomlist, ele[2], ele[3] ]);
                        print '\t', ele[0], '\t', ele[1], '\t', ele[2], '\t', ele[3]; 
                elif name == 'angles':
                    for ele in elements:
                        #12    15    17      2    127.0    20.0
                        angles.append([ ele[0] + ndx_atomlist, ele[1] + ndx_atomlist, ele[2] + ndx_atomlist, ele[3], ele[4], ele[5] ]);
                        print '\t', ele[0]+ndx_atomlist, '\t', ele[1] +ndx_atomlist, '\t', ele[2] +ndx_atomlist, '\t', ele[3], '\t', ele[4], '\t', ele[5];
                elif name == 'dihedrals':
                    for ele in elements:
                        # 1 2 3 1 2 200.0 200.0
                        dihedrals.append([ ele[0]+ndx_atomlist, ele[1]+ndx_atomlist, ele[2]+ndx_atomlist, ele[3], ele[4], ele[5], ele[6] ])
                        print '\t', ele[0]+ndx_atomlist, '\t', ele[1]+ndx_atomlist, '\t', ele[2]+ndx_atomlist, '\t', ele[3], '\t', ele[4], '\t', ele[5], '\t', ele[6]

            
            print "------------------PRINT FINAL SYSTEM INFO-------------------------------"
            print atoms
            print bonds
            print constraints
            print angles
            print dihedrals
            print "------------------PRINT FINAL SYSTEM INFO END ---------------------------"

    
        ####################
        ####################
        ####################
        self.DnpAtomList=mol
        self.DnpAtomITPList=atoms
        self.DnpBondList=bonds
        self.DnpAngleList=angles
        self.DnpConstraintsList=constraints
        self.DnpdihedralList=dihedrals


        print "========================="
        print "=SYSTEM GENERATE SUCCESS="
        print "========================="

    def DNPwritePDB(self):
        #PRINT PDB
        filename="DNP.pdb"
        try:
            fp = open(filename , 'w')
        except:
            print "Error: No such file: "+filename
            exit(1)
        for j,p in enumerate(self.DnpAtomList):
            atom = PDB_module.Atom_class(AtomSerial=p[0], AtomName=p[2], ResidueName=p[3], ResidueSerial=p[1], AtomCoorX=p[4], \
                                         AtomCoorY=p[5], AtomCoorZ=p[6])	
            fp.write("%s\n" %(atom.atom_2_PDBformat()))
        fp.close()

        print "========================="
        print "=    PDB WRITE ENDED    ="
        print "========================="


    def DNPwriteITP(self):
        #PRINIT ITP
        filename="DNP.itp"
        try:
            fp = open(filename , 'w')
        except:
            print "Error: No such file: "+filename
            exit(1)

        #tile
        fp.write("; MARTINI (martini22) Coarse Grained topology file for NP with grafts\n")
        fp.write("\n")
        fp.write("[ moleculetype ]\n")
        fp.write("; name nrexcl\n")
        fp.write("DNP     3\n")
        fp.write("\n")

        line=''

        #ATOM part
        fp.write("[ atoms ]\n")
        fp.write("; id    type    resnr   residu  atom    cgnr    charge\n")
        for j,p in enumerate(self.DnpAtomITPList):
            line = GRO_TOP.Atom2ITP(nr=p[0], atomtype=p[1], resnr=p[2], residu=p[3], atomlabel=p[4], cgnr=p[5], charge=p[6])
            fp.write("%s\n" %(line))
        fp.write("\n")

        #BOND part
        fp.write("[ bonds ]\n")
        fp.write(";  i  j         funct   length  force.c.\n")
        for j,p in enumerate(self.DnpBondList):
            line = GRO_TOP.Bond2ITP(ai=p[0], aj=p[1], funct=p[2], c0=p[3], c1=p[4])
            fp.write("%s\n" %(line))
        fp.write("\n")

        #CONSTRAINTS part
        fp.write("[ constraints ]\n")
        fp.write(";  i  j  k      funct   length\n")
        for j,p in enumerate(self.DnpConstraintsList):
            line = GRO_TOP.Constraints2ITP(ai=p[0], aj=p[1], funct=p[2], c0=p[3])
            fp.write("%s\n" %(line))
        fp.write("\n")

        #ANGLE part
        fp.write("[ angles ]\n")
        fp.write(";  i  j  k      funct   angle   force.c.\n")
        for j,p in enumerate(self.DnpAngleList):
            line = GRO_TOP.Angle2ITP(ai=p[0], aj=p[1], ak=p[2], funct=p[3], c0=p[4], c1=p[5])
            fp.write("%s\n" %(line))
        fp.write("\n")

        #Dihdral part
        fp.write("[ dihdrals ]\n")
        fp.write(";  i  j  k  l   funct   angle   force.c.\n")
        for j,p in enumerate(self.DnpdihedralList):
            line = GRO_TOP.Dihedral2ITP(ai=p[0], aj=p[1], ak=p[2], al=p[3], funct=p[4], c0=p[5], c1=p[6], c2=p[7])
            fp.write("%s\n" %(line))
        fp.write("\n")


        fp.close()

        print "========================="
        print "=    ITP WRITE ENDED    ="
        print "========================="



def generateSystem(Ncore=40, Ngraft=2, Bndgraft=1.0, Radcore=8.0, pdbfile="pdb.pdb", refatom1=1, refatom2=20, itpfile="itp.itp") : 
    sdnp = DnaNanoParticle(NpBall = Ncore, Ndna=Ngraft, RadiurBall = Radcore, BondSP = Bndgraft, \
    	pdbname = pdbfile, fatom1 = refatom1, fatom2 = refatom2, itpname = itpfile)
    sdnp.DNPwritePDB()
    sdnp.DNPwriteITP()
   


if __name__=="__main__":
    ###################################
    # Lsstrand: The length of single strand DNA
    # Lscpaer : The length of double strand DNA
    # Llinker : The length of linkers
    # Bndsstrand: bond length betwee two single strand DNA beads
    # Bndgraft  : bond length betwee two double strand DNA beads
    # BndFLs    : bond length betwee two linkder beads
    ###################################
    generateSystem(Ncore=250, Ngraft=1, Bndgraft=3.0, Radcore=15.0, pdbfile="lys10_val6_md.pdb", refatom1=0, refatom2=40, itpfile="lys10_val6.b.itp")
 
