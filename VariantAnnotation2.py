#Developed modules
import va
import mvi
import vep
import hpa
import uniprot
#
def generatefilterHGVS():
    #Variables
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
                listaffected[li].append(idHGVS.strip('\n'))
        #For control:
        if cat == 'c':
            listcontrol.append([])
            li = len(listcontrol)- 1
            for line in inputfile:
                listcontrol[li].append(idHGVS.strip('\n'))
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
    outputHGVS(listfiltered, "candidate")
    outputHGVS(lglvariant, "GLvariant")

#Variables
listanno = []
listmvi = []
listvep = []
listhpa = []
listuni = []

#Main Menu
exit_program = False
while exit_program == False:
    #Display menu
    print('1) Import data')
    print('2) Export data')
    print('3) Filter')
    print('4) Other')
    print('5) Exit')
    #Accept input
    option = input('Select option: ')
    #Input check
    while (option != '1'
           and option != '2'
           and option != '3'
           and option != '4'
           and option != '5'
           ):
        option = input('Invalid selection; please select one of the options: ')
    print('')
    #1) Import data
    if option == '1':
        print('Yes')
    #2) Export data
    elif option == '2':
        print('Yes')
    #3) Filter
    elif option == '3':
        print('1) Filter affected and control')
        print('2) Filter by consequence terms')
        print('3) Filter by allele frequency')
        print('4) Filter by brain expression')
        print('5) Exit')
        option = input('Select option: ')
        #Input check
        while (option != '1'
               and option != '2'
               and option != '3'
               and option != '4'
               and option != '5'
               ):
            option = input('Invalid selection; please select one of the options: ')
        print('')
        if option == '1':
            generatefilterHGVS()
            print('Filtering step complete.')
        elif option == '2':
            print('Option2')
        elif option == '3':
            print('Option3')
        elif option == '4':
            print('Option4')
        elif option == '5':
            print('Returning to main menu....')
    #4) Other
    elif option == '4':
        print('1) Pull HGVS IDs from VCF files')
        print('2) Get sequences from Ensembl')
        print('3) Annotate using data source')
        print('4) Exit')
        option = input('Select option: ')
        #Input check
        while (option != '1'
               and option != '2'
               and option != '3'
               and option != '4'
               ):
            option = input('Invalid selection; please select one of the options: ')
        print('')
    #5) Exit
    elif option == '5':
        check = input('Are you sure you want to exit? (y/n)')
        if check == 'y':
            exit_program = True
        else:
            print('Returning to menu' + '\n')
