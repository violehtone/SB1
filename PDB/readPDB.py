#!/usr/bin/python

"""
Structural Bioinformatics assignment 1 - Phi-psi angles

When you finish the assignment, the script should create an
output file containing the phi and psi angles and secondary
structure assignment of each residue. Please ONLY modify the
code in the three indicated blocks and do NOT use additional
python packages. Use spaces instead of tabs. Do not round
down in any calculation step.

To run, make 'PDB' your working directory and use:

> python3 readPDB.py pdb_filename.txt
"""

# Packages
from sys import argv
import os
from math import sqrt, atan2, degrees


# Vector functions that we need to calculate the angles
def vectorSubtract(v1, v2):
    return [a - b for a, b in zip(v2, v1)]


def unitVector(v1):
    return [x / magnitude(v1) for x in v1]


def dot_product(v1, v2):
    """ Calculate the dot product of two vectors """
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]


# print(dot_product([1, 2, 3], [1, 3, 2]))


def cross_product(v1, v2):
    """ Calculate the cross product of two vectors """
    i = v1[1] * v2[2] - v1[2] * v2[1]
    j = v1[2] * v2[0] - v1[0] * v2[2]
    k = v1[0] * v2[1] - v1[1] * v2[0]
    return [i, j, k]


# print(cross_product([1, 2, 3], [1, 3, 2]))


def magnitude(v):
    """ Calculate the size of a vector """
    return sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)


# print(magnitude([1, 2, 2]))

# PDB file parser


def readPDB(PDB_file):
    """ Reads a PDB file and stores the atom
    coordinates and amino acid types of the protein """
    # open the file
    try:
        f = open(PDB_file, 'r')
    except:
        print('Error: cannot open', PDB_file)

    # dictionaries to store the output
    # pdb atom coordinates:
    #     pdbcoord[chain][residue_number][atom_type] = coordinates
    pdbcoord = {}
    # residue type per chain and residue number (i.e. store the sequence)
    #     pdbseq[chain][resnum] = restype
    pdbseq = {}

    # parse each line in the file
    for line in f:
        # remove whitespace at the end of the line
        line = line.strip()
        # only parse the lines containing atom coordinates
        if line[:4] == 'ATOM':
            # ATOM type (e.g. C-alpha)
            atom_type = line[12:16].strip()
            # AMINO ACID type (e.g. alanine)
            aa_type = line[17:20].strip()
            # residue number
            res_num = int(line[22:26])
            # Protein chain
            chain = line[21]
            # coordinates
            xcoord = float(line[30:38])
            ycoord = float(line[38:46])
            zcoord = float(line[46:54])

            # if chain does not exists create new entry
            if not chain in pdbcoord:
                pdbcoord[chain] = {}
                pdbseq[chain] = {}
            # if resnum does not exists create new entry
            if not res_num in pdbcoord[chain]:
                pdbcoord[chain][res_num] = {}

            # store coordinates as a vector
            pdbcoord[chain][res_num][atom_type] = [xcoord, ycoord, zcoord]
            # store sequence
            pdbseq[chain][res_num] = aa_type

    # close file
    f.close()

    # return dictionaries
    return pdbcoord, pdbseq


# print(readPDB('Data\\1TIM.pdb'))

# THE FOLLOWING THREE FUNCTIONS ARE THE ONES YOU NEED
# TO EDIT FOR THE ASSIGNMENT. ONLY EDIT THE INDICATED
# BLOCKS


def calculateDihedral(a1, a2, a3, a4):
    """ Calculates the normal vector of the planes
    defined by four atom coordinates """
    # START CODING HERE
    # calculate normal vectors to the planes defined by a1,a2,a3 and a2,a3,a4
    # you may use the functions "cross_product","dot_product" and "magnitude" defined above
    # you can also use the python math function "atan2" and "degrees"

    #      vector             phi          psi
    #   -> v1 = a1 -> a2    (C->N2)   / (N2->Ca2)
    #   -> v2 = a2 -> a3    (N2->Ca2) / (Ca2->C2)
    #   -> v3 = a2 -> a3    (N2->Ca2) / (Ca2->C2)
    #   -> v4 = a3 -> a4    (Ca2->C2) / (C2->N3)

    #                         a1   a2    a3    a4
    # phi = calculateDihedral(c,   n2,   ca2,  c2)
    # psi = calculateDihedral(n2,  ca2,  c2,   n3)

    #define any two vectors of the 3 points that define a plane
    v1 = vectorSubtract(a1, a2)
    v2 = vectorSubtract(a2, a3)
    v3 = vectorSubtract(a2, a3) #a3,a2
    v4 = vectorSubtract(a3, a4)

    #v2=v3, u(v2)


    # define normals to the planes
    n1 = cross_product(v1, v2)
    n2 = cross_product(v3, v4)

    # Calculate angle between the normals
    cos = dot_product(n1, n2) / (magnitude(n1) * magnitude(n2))
    u = unitVector(v2)
    #sin = magnitude(cross_product(n1,n2)) / (magnitude(n1) * magnitude(n2))
    sin = dot_product(cross_product(n1,n2), u) / (magnitude(n1) * magnitude(n2))
    dihedral = degrees(atan2(sin, cos))

    # END CODING HERE
    return dihedral


# print(caculateDihedral([1, 9, 2], [3, 2, 1], [2, 4, 7], [8, 2, 5]))


def assign_ss(phi, psi):
    """ Assign a secondary structure type based on the phi
    and psi angles of a residue """
    # START CODING HERE
    # for code checking purposes use the terms "loop", "alpha" or "beta"

    secondary_structure = ""

    if psi > 0 and phi > 0:
        secondary_structure = "alpha"  # left handed
    elif psi > 0 and phi < 0:
        secondary_structure = "beta"
    elif psi < 0 and phi < 0:
        secondary_structure = "alpha"  # right handed
    else:
        secondary_structure = "loop"

    # END CODING HERE
    return secondary_structure


# print(assign_ss(60, 25))


def print_phi_psi(pdbcoord, pdbseq, outfile):
    """ given the PDB coordinates, calculate the dihedral
    angles of all the residues, assign secondary structure
    types and write them into an output file """
    f = open(outfile, 'w')

    # get the chains from the PDB file
    list_chains = sorted(pdbcoord.keys())

    for chain in list_chains:
        # get the sorted residue numbers from the pdbcoord dictionary
        list_residue_numbers = sorted(pdbcoord[chain].keys())
        for res_num in list_residue_numbers:
            # if certain residues are missing in the PDB file, you will
            # get a KeyError. Make sure your program does not crash, but
            # gives a warning when this happens
            try:
                # START CODING HERE
                #skip first and last residue
                if(res_num == list_residue_numbers[0] or res_num == list_residue_numbers[-1]):
                    continue

                res_num0 = list_residue_numbers[list_residue_numbers.index(res_num)-1]
                res_num2 = list_residue_numbers[list_residue_numbers.index(res_num)+1]

                # define the atoms' coordinates used in phi/psi
                # pdbcoord[chain][res_num][atom_type] => [x,y,z]
                c0 = pdbcoord[chain][res_num0]['C']  # previous residue C
                n = pdbcoord[chain][res_num]['N']  # current residue N
                ca = pdbcoord[chain][res_num]['CA']  # current residue Ca
                c = pdbcoord[chain][res_num]['C']  # current residue C
                n2 = pdbcoord[chain][res_num2]['N']  # next residue N

                phi = calculateDihedral(c0,n, ca, c)
                psi = calculateDihedral(n, ca, c, n2)
                ss = assign_ss(phi,psi)

            # END CODING HERE
            except KeyError:
                print('WARNING: KeyError:', KeyError,
                      'in residue', chain, res_num)

            # get amino acid
            aa_type = pdbseq[chain][res_num]
            # write into output file
            print(chain, res_num, aa_type, phi, psi, ss, file=f)
    f.close()
    print('written:', outfile)


def main():
    # input PDB file
    f_in = argv[1]
    f_out = 'Output/phi_psi.txt'
    # read PDB file
    pdbcoord, pdbseq = readPDB(f_in)
    print_phi_psi(pdbcoord, pdbseq, f_out)
    # for testing
    # for i in ['1TIM', '3PG8']:
    #     f_in = 'student/{}.pdb'.format(i)
    #     print(f_in)
    #     f_out = 'student/output/phi_psi_{}.txt'.format(i)
    #
    #     # read PDB file
    #     pdbcoord, pdbseq = readPDB(f_in)
    #     print_phi_psi(pdbcoord, pdbseq, f_out)


if __name__ == '__main__':
    main()
