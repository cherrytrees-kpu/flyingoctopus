"""
VariantAnnotation.py - handles variant annotation tasks

Author: Michael Ke

"""
import json
import sys
import requests
import myvariant
import time
import va
import datetime

def fullroutine():
    #Get unique job ID name
    now = datetime.datetime.now()
    timelabel = now.strftime('%y%m%d-%H%M%S')

    #Import candidate variants
    lCandidate = va.importfHGVS()

    #MVI Annotation
    startmvi = time.time()
    lmvi = va.annotmvi(lCandidate)
    endmvi = time.time()
    print('Total MVI annotation run time: ' + str(endmvi - startmvi) + '\n')
    #Export MVI
    mvifname = 'mvi_annotation_' + timelabel + '.txt'
    va.exportanno(lmvi, mvifname)
    #VEP annotation
    startvep = time.time()
    lvep = va.annotvep(lCandidate)
    endvep = time.time()
    print ('Total VEP annotation run time: '+  str(endvep - startvep) + '\n')
    #Export VEP
    vepfname = 'vep_annotation_' + timelabel + '.txt'
    va.exportanno(lvep, vepfname)
    #Export combined annotations
    lanno = va.combineanno(lmvi, lvep)
    annofname = 'annotated_mutations_' + timelabel + '.txt'
    va.writeanno(lanno, annofname)
    #Filtervariant
    va.filtervariant(lanno, timelabel)
    print('Complete')
##### Program Start ##################################################################
#GLOBAL VARIABLES
EXIT_PROGRAM = False
#LIST_ANNO - holds all annotations
LIST_ANNO = []
#LIST_VEP - holds VEP annotations
LIST_VEP = []
#LIST_MVI - holds myvariant.info annotations
LIST_MVI = []

while EXIT_PROGRAM == False:
    print('1) Import HGVS from .vcf files')
    print('2) Generate filtered HGVS file')
    print('3) Annotation')
    print('4) Import data')
    print('5) Export data')
    print('6) Process data')
    print('7) Full Routine')
    print('8) Exit (modules still loaded)')

    option = int(input('Select option: '))

    while (option != 1
           and option != 2
           and option != 3
           and option != 4
           and option != 5
           and option != 6
           and option != 7
           and option != 8
           ):
        option = int(input('Invalid selection; please select one of the options: '))

    print('')
    #1) Import HGVS from .vcf files
    if option == 1:
        try:
            filename = input('Please enter the .vcf file to extract HGVS from: ')
        except IoError:
            print ('File not found.')
            filename = input('Re-enter filename: ')
        va.vcftoHGVS(filename)

    #2) Generate filtered HGVS file
    if option == 2:
        #Variables
        start = time.time()
        listaffected = []
        listcontrol = []
        numfiles = int(input('How many files individuals would you like to use: '))
        i = 0

        #Import the number of files that will be analyzed
        while i < numfiles:
            #Open the HGVS ID file
            try:
                filename = input('Please enter the name of the file listing HGVS IDs to be imported: ')
                inputfile = open(filename, 'r')
            except IOError:
                print('File not found.')
                filename = input('Re-enter filename: ')

            #Ask user if the HGVS ID file is from an affected or control
            cat = input('Enter "a" for HGVS IDs from affected, "c" from controls: ')
            while cat != 'a' and cat != 'c':
                cat = input('Please enter "a" for affected, and "c" for control: ')

            #For affected:
            if cat == 'a':
                listaffected.append([])
                li = len(listaffected) - 1
                for idHGVS in inputfile:
                    listaffected[li].append(line.strip('\n'))
            #For control:
            if cat == 'c':
                listcontrol.append([])
                li = len(listcontrol)- 1
                for line in inputfile:
                    listcontrol[li].append(line.strip('\n'))
            print(filename + ' successfully imported.')
            i = i + 1

        #Perform filtering
        listfiltered= va.filteraffected(listaffected, listcontrol)
        list_gl = []

        #Remove any GL variants from the filtered list
        j = 0
        while j < len(listfiltered):
            if 'chrGL' in listfiltered[j]:
                list_gl.append(listfiltered[j])
                del listfiltered[j]
                j = j - 1
            j = j + 1

        va.outputHGVS(lCandidate, "candidate")
        va.outputHGVS(lglvariant, "GLvariant")
        end = time.time()
        print('Total run time: ' + str(end - start) + '\n')

    elif option == 3:
        exitannotation = False
        while exitannotation == False:
            print ('1) myvariant.info annotation')
            print ('2) VEP annotation')
            print ('3) Return to main menu')

            optionannotation = int(input('Select option: '))

            while (optionannotation != 1
                   and optionannotation != 2
                   and optionannotation != 3):
                optionannotation = int(input('Invalid selection; please select one of the options: '))

            if optionannotation == 1:
                #Open HGVS ID file that will be annotated
                try:
                    filename = input('Please enter the name of the file listing HGVS IDs to be imported: ')
                    inputfile = open(filename, 'r')
                except IOError:
                    print('File not found.')
                    filename = input('Re-enter filename: ')
                listHGVS = va.importHGVS(filename)
                #Perform MVI annotation
                start = time.time()
                LIST_MVI = va.annotmvi(listHGVS)
                end = time.time()
                print('Total run time: ' + str(end - start) + '\n')

            elif optionannotation == 2:
                #Open HGVS ID file that will be annotated
                try:
                    filename = input('Please enter the name of the file listing HGVS IDs to be imported: ')
                    inputfile = open(filename, 'r')
                except IOError:
                    print('File not found.')
                    filename = input('Re-enter filename: ')
                listHGVS = va.importHGVS(filename)
                #Perform VEP annotation
                start = time.time()
                LIST_VEP = va.annotvep(listHGVS)
                end = time.time()
                print('Total run time: ' + str(end - start) + '\n')

            elif optionannotation == 3:
                exitannotation = True
                print('')

    elif option == 4:
        exitimport = False
        while exitimport == False:
            print ('1) Import myvariant.info annotation data')
            print ('2) Import VEP annotation data')
            print ('3) Import program-generated data')
            print ('4) Return to main menu')

            optionimport = int(input('Select option: '))

            while (optionimport != 1
                   and optionimport !=2
                   and optionimport !=3
                   and optionimport !=4):
                optionimport = int(input('Invalid selection; please select one of the options: '))

            if optionimport == 1:
                LIST_MVI = va.importanno()
                print('Data imported.' + '\n')
            elif optionimport == 2:
                LIST_VEP = va.importanno()
                print('Data imported.' + '\n')
            elif optionimport == 3:
                LIST_ANNO = va.importanno()
                print('Data imported.' + '\n')
            elif optionimport == 4:
                exitimport = True
                print('')

    elif option == 5:
        exitexport = False
        while exitexport == False:
            print('1) Export myvariant.info annotation data')
            print('2) Export VEP annotation data')
            print('3) Export tab-delimited data')
            print('4) Return to main menu')

            optionexport = int(input('Select option: '))

            while (optionexport != 1
                   and optionexport != 2
                   and optionexport != 3
                   and optionexport != 4):
                optionexport = int(input('Invalid selection; please select one of the options: '))

            if optionexport == 1:
                va.exportanno(LIST_MVI, 'myvariantannotation.txt')
                print('Export completed.' + '\n')
            elif optionexport == 2:
                va.exportanno(LIST_VEP, 'vepannotation.txt')
                print('Export completed.' + '\n')
            elif optionexport == 3:
                LIST_ANNO = va.combineanno(LIST_MVI, LIST_VEP)
                va.export(LIST_ANNO, 'annotated_mutations.txt')
            elif optionexport == 4:
                exitexport = True
                print('')

    elif option == 6:
        listfiltered = va.filtervariant(LIST_ANNO)

    elif option == 7:
        fullroutine();

    elif option == 8:
        check = input('Are you sure you want to exit? (y/n)')
        if check == 'y':
            EXIT_PROGRAM = True
        else:
            print('Returning to menu' + '\n')

print ("Terminating program...")
