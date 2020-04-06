#Developed modules
import mvi
import vep
import hpa
import uniprot

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
def main():
    print('Variant Annotation Program V1')
    #Lists of variants
    list_irrelevant = []
    list_highfreq = []
    list_notexpressedbrain = []
    list_candidate = []
    #List of data
    list_mvi = []
    list_vep = []
    list_hpa = []
    list_uni = []

    #Get HGVS IDs
    filename = input('Please enter name of file containing HGVS IDs of variants to analyze: ')
    listhgvs = importHGVS(filename)
    #Begin annotation and filtering
    for id in listhgvs:
        filtered = False
        while filtered == False:
            data_variant = {'_id':id}
            #Get VEP annotations
            vepanno = vep.getvep(id)
            #Parse VEP annotations into data structure
            if vepanno is not None:
                vep.parsevepanno(data_variant, vepanno)
            mvianno = mvi.getmvi(id)


if __name__ == "__main__":
    main()
