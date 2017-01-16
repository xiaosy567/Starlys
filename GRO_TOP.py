
#/usr/bin/env python
import sys 
import os
import string
import re

'''
read pdb
'''
def Read_pdb(filename):
    '''
    Read in pdb file and output the atom list.
    '''
    try:
        fp=open(filename,'r')
    except IOError, e:
        print e
        sys.exit()
    all_lines=fp.readlines()
    
    atom=[]
    for line in all_lines:
        if "ATOM" in line:
            eles = line.split();
            atom.append([eles[0], int(eles[1]), eles[2], eles[3], int(eles[5]), float(eles[6]), float(eles[7]), float(eles[8]), float(eles[9]), float(eles[10]) ])

    print "==ouput PDB fire for checking=="
    for j,p in enumerate(atom):
        for i, ln in enumerate(p):
            print ln, "\t",
        print
    print "==ouput PDB fire for checking=="
    
    # print atom
    return atom


def Read_itp(filename):
    '''
    Read in pdb file and output the atom list.
    '''

    print "==READIND ITP FILE=="

    try:
        fp=open(filename,'r')
    except IOError, e:
        print e
        sys.exit()

    moleculetype=[]
    atoms=[]
    bonds=[]
    angles=[]
    constraints=[]
    dihedrals=[]

    all_lines=fp.readlines()
    lines=[];
    for line in all_lines:
        if line[0:1] != ";" and len(line)>1:
            lines.append(line[:-1]);
            print line[:-1];

    pointer=0;
    dict ={'moleculetype':0, 'atoms':1, 'bonds':2, 'angles':3, 'constraints':4, 'dihedrals':5}
    print "======="
    for j, p in enumerate(lines):
        if p[0:1] =="[":
            element=p.split()[1];
            pointer=dict[element];
            # print pointer;
            continue;

        elif pointer == 0: # moleculetype
            eles=p.split();
            try:
                eles.remove(';')
            except:
                pass
            moleculetype.append(eles)

        elif pointer == 1: # atoms
            eles=p.split();
            try:
                eles.remove(';')
            except:
                pass
            #format: 1    Qd     1   LYS    BB     1  1.0000
            atoms.append([ int(eles[0]), eles[1], int(eles[2]), eles[3], eles[4], int(eles[5]), float(eles[6]) ]);

        elif pointer == 2: #bonds;
            eles=p.split();
            try:
                eles.remove(';')
            except:
                pass
            #format: 1     4      1   0.35000  1250
            bonds.append([ int(eles[0]), int(eles[1]), int(eles[2]), float(eles[3]), float(eles[4]) ]);

        elif pointer == 3: #angles;
            eles=p.split();
            try:
                eles.remove(';')
            except:
                pass
            #format: 1     4     7      2    127.0    20.0
            angles.append([ int(eles[0]), int(eles[1]), int(eles[2]), int(eles[3]), float(eles[4]), float(eles[5]) ]);

        elif pointer == 4: #constraints;
            eles=p.split();
            try:
                eles.remove(';')
            except:
                pass
            # format: 15    16      1   0.26500
            constraints.append([ int(eles[0]), int(eles[1]), int(eles[2]), float(eles[3]) ]);

        elif pointer == 5: #dihedrals;
            eles=p.split();
            try:
                eles.remove(';')
            except:
                pass
            dihedrals.append([ int(eles[0]), int(eles[1]), int(eles[2]), int(eles[3]), int(eles[4]), float(eles[5]), float(eles[6]) ]);

        else:
            print "!!!XTRA element exists in ITP file, CHECK IT!"
            print "!!!NOT atoms, bonds, angles, angles, dihedrals, constraints!"

    itp_list=[]
    itp_list.append(["moleculetype", moleculetype])
    itp_list.append(["atoms", atoms])
    itp_list.append(["bonds", bonds])
    itp_list.append(["constraints", constraints])
    itp_list.append(["angles", angles])
    itp_list.append(["dihedrals", dihedrals])

    print "==END READING ITP FILE=="

    return itp_list

def print_itp(itplist):
    '''
    print itp file
    '''
    itp_list=itplist
    print itp_list
    for j, p in enumerate(itp_list):
        name=p[0];
        elements=p[1];
        print '[ ', name, ' ]';
        if name == 'moleculetype':
            for ele in elements:
                print '\t', ele[0], '\t', ele[1];
            print ""

        elif name == 'atoms':
            for ele in elements:
                #1    Qd     1   LYS    BB     1  1.0000
                print '\t', ele[0], '\t', ele[1], '\t', ele[2], '\t', ele[3], '\t', ele[4], '\t', ele[5], '\t', ele[6]
            print ""

        elif name == 'bonds':
            for ele in elements:
                #1     4      1   0.35000  1250
                print '\t', ele[0], '\t', ele[1], '\t', ele[2], '\t', ele[3], '\t', ele[4]
            print ""

        elif name == 'constraints':
            for ele in elements:
                #15    16      1   0.26500
                print '\t', ele[0], '\t', ele[1], '\t', ele[2], '\t', ele[3]
            print ""

        elif name == 'angles':
            for ele in elements:
                #12    15    17      2    127.0    20.0
                print '\t', ele[0], '\t', ele[1], '\t', ele[2], '\t', ele[3], '\t', ele[4], '\t', ele[5]
            print ""

        elif name == 'dihedrals':
            for ele in elements:
                # 1 2 3 1 2 200.0 200.0
                print '\t', ele[0], '\t', ele[1], '\t', ele[2], '\t', ele[3], '\t', ele[4], '\t', ele[5], '\t', ele[6]
            print ""

def Atom2ITP(nr=1, atomtype='C', resnr=1, residu='UREA', atomlabel='C1', cgnr=1, charge=0.683):
    '''
    TEMPLATE:
       1      C          1      UREA     C1   0.683
    OUT E.G.:
       1      C          1      UREA     C1   0.683
    '''
    s="%6d%8s%8d%8s%8s%8d%8.3f" %(nr, atomtype, resnr, residu, atomlabel, cgnr, charge)  
    return s

def Bond2ITP(ai=3, aj=4, funct=1, c0=0.1, c1=37500):
    '''
    TEMPLATE:
        1     4      1   0.35000  1250
    OUT.E.G.
        1     4      1   0.35000  1250
    '''
    s="%6d%6d%6d%10.3f%10.3f" %(ai, aj, funct, c0, c1)
    return s

def Constraints2ITP(ai=3, aj=4, funct=1, c0=0.1):
    '''
    TEMPLATE
        7     8      1   0.26500
    '''
    s="%6d%6d%6d%10.3f" %(ai, aj, funct, c0)
    return s

def Angle2ITP(ai=1, aj=3, ak=4, funct=1, c0=120.0, c1=290.0):
    '''
    TEMPLATE
        1  3 4 1 120.0  292.8
    '''
    s="%6d%6d%6d%6d%10.3f%10.3f" %(ai, aj, ak, funct, c0, c1)
    return s

def Dihedral2ITP(ai=1, aj=3, ak=4, al=5, funct=1, c0=120.0, c1=29.0, c2=2.0):
    '''
    TEMPLATE
        1  3 4  5  1 120.0  29.2, 2.0
    '''
    s="%6d%6d%6d%6d%6d%10.3f%10.3f" %(ai, aj, ak, al, funct, c0, c1, c2)
    return s

# if __name__=="__main__":
#     # filename_pdb="lys10_val6_md.pdb"
#     # atom_list=Read_pdb(filename_itp)
#     filename_itp="lys10_val6.b.itp"
#     itp_list=Read_itp(filename_itp)
#     print_itp(itp_list)