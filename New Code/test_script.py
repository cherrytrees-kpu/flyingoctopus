import myvariant
import requests
import json

def importHGVS(filename):
    
    """
    importHGVS - imports HGVS into the program
    Parameters: none
    Return: list containing the imported HGVS ids
    """

    inputfile_open = False
    listHGVS = []

    #Open the file
    while inputfile_open == False:
        try:
            inputfile = open(filename, 'r')
            inputfile_open = True
        except IOError:
            print('File not found.\n')

            filename = input('Re-enter filename: ')

    #Import data into listHGVS
    for line in inputfile:
        listHGVS.append(line.strip('\n'))

    inputfile.close()

    return listHGVS

def exportanno(listanno, filename):
    """
    exportanno - exports raw annotation data
    Parameters:
    listanno - list - list of dictionaries containing raw JSON data
    filename - string - name of the file to output to
    Returns: none; outputs to file with filename
    """
    #Create output file
    outputfile = open(filename, 'w')

    #For each annotation in the annotation list, dump the JSON data to file
    for anno in listanno:
        outputfile.write(json.dumps(anno))
        #Check to see if the last element has been reached before entering a
        #newline character
        if listanno.index(anno) != (len(listanno)-1):
            outputfile.write('\n')
    outputfile.close()
    print('Export completed.' + '\n')

def annotateMVI(HGVS):

    """
    annotmvi - accepts HGVS, returns data on mutation
    Parameters: listHGVS: list of HGVS IDs to retrieve annotations for
    Return: dictionary/list containing the json data from myvariant.info. Returns none if no entry is found.
    """
    mv = myvariant.MyVariantInfo()

    #For each HGVS ID in the list, retrieve the annotation
    annoMVI = mv.getvariant(HGVS, fields = [
        'cadd', 
        'clinvar',
        'dbnsfp',
        'dbsnp',
        'gnomad_exome',
        'gnomad_genome',
        ])

    return annoMVI

if __name__ == "__main__":
    filename = input('Please enter HGVS filename: ')
    listHGVS = importHGVS(filename)
    
    listMVI = []

    for HGVS in listHGVS: 
        listMVI.append(annotateMVI(HGVS))
        if listHGVS.index(HGVS)%1000 == 0: 
            print (str(listHGVS.index(HGVS)) + 'out of ' + str(len(listHGVS)) + 'completed')

    exportanno(listMVI, 'mvi_annotation.txt')