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
import PDB_module
import GRO_TOP

def check_BOND(list1=[" a ", " b"], list2=["b", "c "]):
	'''
	same bond return 1
	otherwise return 0
	'''
	
	pointer = 0

	list1_st=[]
	list2_st=[]
	for ele in list1:
		list1_st.append(string.strip(ele))

	for ele in list2:
		list2_st.append(string.strip(ele))

	if ((list1_st[0]==list2_st[0]) and (list1_st[1] == list2_st[1])) :
		pointer=1
	elif ((list1_st[0]==list2_st[1]) and (list1_st[1] == list2_st[0])):
			pointer = 1

	return pointer


def check_ANGLE(list1=[" a ", " b", " c "], list2=["b", "c ", "d"]):
	'''
	same agnle return 1
	otherwise return 0
	'''

	pointer = 0
	list1_st=[]
	list2_st=[]
	for ele in list1:
		list1_st.append(string.strip(ele))

	for ele in list2:
		list2_st.append(string.strip(ele))
	
	if ((list1_st[0]==list2_st[0]) and (list1_st[1] == list2_st[1]) and (list1_st[2] == list2_st[2])) :
		pointer=1
	elif ((list1_st[0]==list2_st[2]) and (list1_st[1] == list2_st[1]) and (list1_st[2] == list2_st[0])):
			pointer = 1

	return pointer



class itpclass:
	'''
	for itp convert to psf (wataru format PSF)

	NOTE: 1. constraints/dihedrals will be !ignored! during the convert.
			The resulted PSF can only be used for vmd view, combined with my psftemplate program. 
		2.  the bond and angle types are decided based on the atom type (which may not correct in reality)!!!

	'''
	def __init__(self, itpfilename = "lys10_val6.b.itp", psffile = "lys10_val6.b.psf"):

		self.atomsLIST=[]
		self.bondsLIST=[]
		self.anglesLIST=[]
		self.constraintsLIST=[]
		self.dihedralsLIST=[]
		self.bondtypeLIST=[]
		self.angletypeLIST=[]

		print "hello world!"
		print itpfilename

		##read_itp
		toplist=GRO_TOP.Read_itp(itpfilename)
		print "== itp info =="
		print toplist
		print "== itp info =="

		for j, p in enumerate(toplist):
			name=p[0];
			elements=p[1];
			if name == 'moleculetype':
				moleculetype=elements[0]
			elif name == 'atoms':
				atoms = elements
			elif name == 'bonds':
				bonds = elements
			elif name == 'constraints':
				constraints = elements
			elif name == 'angles':
				angles = elements
			elif name == 'dihedrals':
				dihedrals = elements

		BOND_typelist=[]
		ANGLE_typelist=[]

		##BOND.
		for j,p in enumerate(bonds):
			print p;
			b_aindex1=p[0]-1
			b_aindex2=p[1]-1
			# print b_aindex1, b_aindex2
			if j == 0:
				print atoms[b_aindex1][1], atoms[b_aindex2][1]
				BOND_typelist.append([ atoms[b_aindex1][1], atoms[b_aindex2][1] ])
			else:
				tmplist=[ atoms[b_aindex1][1], atoms[b_aindex2][1] ]
				pointer=0
				for k, ele in enumerate(BOND_typelist):
					pointer = check_BOND(tmplist, ele)
					if pointer==1:
						break
					elif (pointer == 0 and k==(len(BOND_typelist)-1) ):
						BOND_typelist.append(tmplist)
						print tmplist

		for j,p in enumerate(angles):
			print p;
			ag_aindex1=p[0]-1
			ag_aindex2=p[1]-1
			ag_aindex3=p[2]-1
			if j == 0:
				print atoms[ag_aindex1][1], atoms[ag_aindex2][1], atoms[ag_aindex3][1]
				ANGLE_typelist.append([ atoms[ag_aindex1][1], atoms[ag_aindex2][1], atoms[ag_aindex3][1] ])
			else:
				tmplist=[ atoms[ag_aindex1][1], atoms[ag_aindex2][1], atoms[ag_aindex3][1] ]
				pointer=0
				for k, ele in enumerate(ANGLE_typelist):
					pointer = check_ANGLE(tmplist, ele)
					if pointer==1:
						break
					elif (pointer == 0 and k==(len(ANGLE_typelist)-1) ):
						ANGLE_typelist.append(tmplist)
						print tmplist

		self.atomsLIST=atoms
		self.bondsLIST=bonds
		self.anglesLIST=angles
		self.constraintsLIST=constraints
		self.dihedralsLIST=dihedrals
		self.bondtypeLIST=BOND_typelist
		self.angletypeLIST=ANGLE_typelist

	def writePSF(self, psffile = "lys10_val6.b.psf"):
		#PRINIT PSF
		filename=psffile
		try:
			fp = open(filename , 'w')
		except:
			print "Error: No such file: "+filename
			exit(1)

		#tile
		fp.write("RESIDNAME  DNP\n")
		fp.write("%s %3d\n" %('NUMATOM', len(self.atomsLIST)))
		fp.write("%s %3d\n" %('NUMBOND', len(self.bondsLIST)))
		fp.write("%s %3d\n" %('NUMANGLE', len(self.anglesLIST)))
		fp.write("%s %3d\n" %('NUMIMPR', 0))

		line=''

		#ATOM part
		for j,p in enumerate(self.atomsLIST):
			line = PDB_module.Atom2PSF(atomHead='ATOM', atomSerial=p[0], atomName=p[4], atomType=p[1], \
								Mass=1.0, Charge=0.0, Unset=0.0)
			fp.write("%s\n" %(line))

		#BOND part
		for j,p in enumerate(self.bondsLIST):
			typeindex=0
			b_aindex1=p[0]-1
			b_aindex2=p[1]-1
			tmplist=[ self.atomsLIST[b_aindex1][1], self.atomsLIST[b_aindex2][1] ]

			pointer=0
			for k, pt in enumerate(self.bondtypeLIST):
				pointer=check_BOND(tmplist, pt)
				if pointer == 1:
					typeindex = k
					break
				elif (pointer == 0 and k==(len(self.bondtypeLIST) - 1) ):
					print "***!!!!*** FIND bond type does not exist! CHECK PROGRAM, please!!!"

			line = PDB_module.Bond2PSF(bondHead='BOND', bondSerial=j+1, bondTypeSerial=typeindex+1, bondIndex1=p[0], bondIndex2=p[1], \
						bondIndex1Type=self.bondtypeLIST[typeindex][0], bondIndex2Type=self.bondtypeLIST[typeindex][1])
			fp.write("%s\n" %(line))

		#ANGLE part 
		for j,p in enumerate(self.anglesLIST):

			typeindex = 0
			ag_aindex1=p[0]-1
			ag_aindex2=p[1]-1
			ag_aindex3=p[2]-1
			tmplist=[ self.atomsLIST[ag_aindex1][1], self.atomsLIST[ag_aindex2][1], self.atomsLIST[ag_aindex3][1] ]

			pointer = 0
			for k, pt in enumerate(self.angletypeLIST):
				pointer = check_ANGLE(tmplist, pt)
				if pointer == 1:
					typeindex = k
				elif (pointer == 0 and k == (len(self.anglesLIST)-1) ):
					print "***!!!!*** FIND angle type does not exist! CHECK PROGRAM, please!!!"

			line = PDB_module.Angle2PSF(angleHead='ANGLE', angleSerial=j+1, angleTypeSerial=typeindex+1, \
						angleIndex1=p[0], angleIndex2=p[1], angleIndex3=p[2], angleIndex1Type=self.angletypeLIST[typeindex][0], \
						angleIndex2Type=self.angletypeLIST[typeindex][1], angleIndex3Type=self.angletypeLIST[typeindex][2])

			fp.write("%s\n" %(line))

		fp.close()

		print "========================="
		print "=    PSF WRITE ENDED    ="
		print "========================="


def generateSystem(itpfilename = "lys10_val6.b.itp", psffile = "lys10_val6.b.psf"):
	myitp=itpclass(itpfilename = "lys10_val6.b.itp")
	myitp.writePSF(psffile = "lys10_val6.b.psf")

if __name__=="__main__":
	# list1=[" a ", " b"]
	# list2=[" b", "a "]
	# a= check_BOND(list1, list2)
	# print a

	# list1=[" a ", " b", "c"]
	# list2=[" c", "b ", "a"]
	# print check_ANGLE(list1, list2)
	generateSystem(itpfilename = "lys10_val6.b.itp", psffile = "lys10_val6.b.psf")

