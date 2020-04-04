#Third-party modules
from consolemenu import *
from consolemenu.items import *
#Developed modules
import va
import mvi
import vep
import hpa
import uniprot

def main():
    def generateimportmenu(main_menu):
        import_menu = ConsoleMenu()
        return SubmenuItem('Import data', import_menu, main_menu)
    def generateexportmenu(main_menu):
        export_menu = ConsoleMenu()
        return SubmenuItem('Export data', export_menu, main_menu)
    def generatefiltermenu(main_menu):
        filter_menu = ConsoleMenu()
        filterstep1 = FunctionItem('Filter by affected and control individuals', filteraffected, [currenthgvs])
        filterstep2 = FunctionItem('Filter by gene consequence', filtercons, [currenthgvs, vepanno])
        filter_menu.append_item(filterstep1)
        filter_menu.append_item(filterstep2)
        return SubmenuItem('Filter variants', filter_menu, main_menu)
    def generateothermenu(main_menu):
        other_menu = ConsoleMenu()
        return SubmenuItem('Other', other_menu, main_menu)
    #Important lists
    list_candidates = []
    currenthgvs = []
    mvianno = []
    vepanno = []
    hpanno = []
    uniprotanno = []
    # Create the menu
    main_menu = ConsoleMenu('VariantAnnotation', 'Program for genetic variant analysis.')
    submenu_import = generateimportmenu(main_menu)
    submenu_export = generateexportmenu(main_menu)
    submenu_filter = generatefiltermenu(main_menu)
    submenu_other = generateothermenu(main_menu)
    main_menu.append_item(submenu_import)
    main_menu.append_item(submenu_export)
    main_menu.append_item(submenu_filter)
    main_menu.append_item(submenu_other)
    main_menu.show()



#Procedural functions
def filteraffected(currenthgvs):
    def gethgvsfiles():
        #VARIABLES
        listaffected = []
        listcontrol = []
        #Get number of files to use in analysis
        try:
            numfiles = int(input('How many individuals would you like to use in analysis: '))
        except ValueError:
            print('Must pass an integer')
        else:
            i = 0
            while i < numfiles:
                try:
                    filename = input('Please enter the name of the file listing HGVS IDs to be imported: ')
                    inputfile = open(filename, 'r')
                except IOError:
                    print('File not found.')
                else:
                    #Ask user if the HGVS ID file is from an affected or control
                    cat = input('Enter "a" for HGVS IDs from affected, "c" from controls: ')
                    while cat != 'a' and cat != 'c':
                        cat = input('Please enter "a" for affected, and "c" for control: ')
                    #Read the file
                    listhgvs = []
                    for idHGVS in inputfile:
                        listhgvs.append(idHGVS.strip('\n'))
                    listhgvs.sort()
                    #For affected:
                    if cat == 'a':
                        listaffected.append(listhgvs)
                    #For control:
                    elif cat == 'c':
                        listcontrol.append(listhgvs)
                    print(filename + ' successfully imported.')
                    i = i + 1
        return listaffected, listcontrol
    listaffected, listcontrol = gethgvsfiles()
    listfiltered = va.filteraffected(listaffected, listcontrol)
    #Remove any GL variants from the filtered list
    j = 0
    while j < len(listfiltered):
        if 'chrGL' in listfiltered[j]:
            list_gl.append(listfiltered[j])
            del listfiltered[j]
            j = j - 1
        j = j + 1
    va.outputHGVS(listfiltered, name = "candidate")
    va.outputHGVS(list_gl, name = "candidate")
    currenthgvs = listfiltered
def filtercons(listhgvs=[], vepanno=[]):
    #Check for HGVS data; if none, import
    if listhgvs == []:
        filename = input ('Please enter name of file containing HGVS IDs of variants to be analyzed: ')
        listhgvs = importHGVS(filename)
    #Check for VEP annotation data


    if vepanno == []:
        print('No data')
    else:
        print('Yes data')

if __name__ == "__main__":
    main()
