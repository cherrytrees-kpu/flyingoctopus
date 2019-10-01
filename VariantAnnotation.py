"""
VariantAnnotation.py v13 - handles variant annotation tasks

Version Date: August 27th, 2019

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

    #Import all of the HGVS files into the Program
    lAffected = []
    lControl = []
    i = 0

    numFiles = int(input('How many files individuals would you like to use: '))
    while i < numFiles:
        va.importHGVS(lAffected, lControl)
        i = i + 1

    lCandidate = va.filteraffected(lAffected, lControl)
    filterfname = "candidate_" + timelabel
    va.outputHGVS(lCandidate, filterfname)

    #MVI Annotation
    startmvi = time.time()
    lmvi = va.annotmvi(lCandidate)
    endmvi = time.time()
    print('Total MVI annotation run time: ' + str(endmvi - startmvi) + '\n')
    #VEP annotation
    startvep = time.time()
    lvep = va.annotvep(lCandidate)
    endvep = time.time()
    print ('Total VEP annotation run time: '+ + str(endvep - startvep) + '\n')
    #Export MVI
    mvifname = 'mvi_annotation_' + timelabel + '.txt'
    va.exportanno(lmvi, mvifname)
    #Export VEP
    vepfname = 'vep_annotation_' + timelabel + '.txt'
    va.exportanno(lvep, vepfname)
    #Export combined annotations
    lanno = va.combineanno(lmvi, lvep)
    annofname = 'annotated_mutations_' + timelabel + '.txt'
    va.writeanno(lanno)
    #Filtervariant
    va.filtervariant(lanno, timelabel)
##### Program Start ##################################################################
#GLOBAL VARIABLES
EXIT_PROGRAM = False
#LIST_ANNO - holds all annotations
LIST_ANNO = []
#LIST_VEP - holds VEP annotations
LIST_VEP = []
#LIST_MVI - holds myvariant.info annotations
LIST_MVI = []
#LIST_PROCESSED - holds filtered variant annotations
LIST_PROCESSED = []

while EXIT_PROGRAM == False:
    print('1) Generate filtered HGVS file')
    print('2) Annotation')
    print('3) Import data')
    print('4) Export data')
    print('5) Process data')
    print('6) Full Routine')
    print('7) Exit (modules still loaded)')

    option = int(input('Select option: '))

    while (option != 1
           and option != 2
           and option != 3
           and option != 4
           and option != 5
           and option != 6
           and option != 7
           ):
        option = int(input('Invalid selection; please select one of the options: '))

    print('')

    if option == 1:
        start = time.time()
        lAffected = []
        lControl = []

        numFiles = int(input('How many files individuals would you like to use: '))

        i = 0
        while i < numFiles:
            va.importHGVS(lAffected, lControl)
            i = i + 1

        lCandidate = va.filteraffected(lAffected, lControl)
        va.outputHGVS(lCandidate)
        end = time.time()
        print('Total run time: ' + str(end - start) + '\n')

    elif option == 2:
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

                listHGVS = va.importfHGVS()

                start = time.time()
                LIST_MVI = va.annotmvi(listHGVS)
                end = time.time()
                print('Total run time: ' + str(end - start) + '\n')

            elif optionannotation == 2:

                listHGVS = va.importfHGVS()
                start = time.time()
                LIST_VEP = va.annotvep(listHGVS)
                end = time.time()
                print('Total run time: ' + str(end - start) + '\n')

            elif optionannotation == 3:
                exitannotation = True
                print('')

    elif option == 3:
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
                LIST_ANNO = va.importmut()
                print('Data imported.' + '\n')
            elif optionimport == 4:
                exitimport = True
                print('')

    elif option == 4:
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
                va.writeanno(LIST_ANNO)
            elif optionexport == 4:
                exitexport = True
                print('')

    elif option == 5:
        listfiltered = va.filtervariant(LIST_ANNO)
        LIST_PROCESSED = va.retrievevep(listfiltered, LIST_VEP)
        va.exportanno(LIST_PROCESSED, 'vepanno_candidates.txt')
        va.transcriptids(LIST_PROCESSED)

    elif option == 6:
        fullroutine();

    elif option == 7:
        check = input('Are you sure you want to exit? (y/n)')
        if check == 'y':
            EXIT_PROGRAM = True
        else:
            print('Returning to menu' + '\n')

print ("Terminating program...")