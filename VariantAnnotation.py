"""
VariantAnnotation.py - handles variant annotation tasks

Author: Michael Ke

"""
#Standard library modules
import json
import sys
import pathlib
import shutil
import requests
import time
import datetime
#Third-party modules
import myvariant
from openpyxl import Workbook
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
#Developed modules
import mvi
import vep
import hpa
import uniprot

#Importing functions
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
def importanno(filename):
    """
    importanno - from a file containing json.dumps annotation data, import it as
    a list of dictionaries
    Parameters: filename - name of file to open
    Return: a list of dictionaries of annotation data
    """

    inputfile_open = False
    listdata = []

    #Open the file
    while inputfile_open == False:
        try:
            inputfile = open(filename, 'r')
            inputfile_open = True
        except IOError:
            print('File not found.\n')
            filename = input('Re-enter filename: ')

    #Import the data
    for line in inputfile:
        listdata.append(json.loads(line))

    inputfile.close()

    return listdata
def vcftoHGVS(filename):
    """
    vcftoHGVS - calls the myvariant.info get_hgvs_from_vcf function
    Parameters: None
    Returns: none; outputs a file containing all of the HGVS IDs
    """
    listHGVS = []

    inputfile_open = False

    while inputfile_open == False:
        try:
            listHGVS = list(myvariant.get_hgvs_from_vcf(filename))
            inputfile_open = True
        except IOError:
            print('File not found.\n')
            filename = input('Re-enter filename: ')

    #Accept name input from user
    name = input('Please enter an identifier for your file: ')
    outputfile = open('HGVS_' + name + '.txt', 'w')

    #Write to file
    print('Writing to file...')
    #for idHGVS in listHGVS:
    length = len(listHGVS)
    i = 0
    while i < (length - 1):
        outputfile.write(listHGVS[i])
        outputfile.write ('\n')
        i = i + 1
    outputfile.write(listHGVS[i])
    print('Done.')

    outputfile.close()
#Exporting functions
def outputHGVS(listHGVS, name = ""):
    """
    outputHGVS - outputs a file with the name 'HGVS_*name.txt' containing the
    list of HGVS ids
    Parameters:
    listHGVS - list - contains list ids
    Returns: nothing; outputs to file 'HGVS_*name.txt'

    """
    i = 0

    #Accept name for output file from user
    if name != "":
        name = input('Please enter an identifier for your file: ')
    outputfile = open('HGVS_' + name + '.txt', 'w')

    #Write to file
    while(len(listHGVS)) > index:
        outputfile.write(listHGVS[i])
        if i != (len(listHGVS)-1):
            outputfile.write('\n')
        i = i + 1

    outputfile.close()
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
def writeanno(listanno, name = "annotated_mutations.txt"):
    """
    writeanno - export annotations to file, tab delimited.
    Parameters:
    listanno - list - contains all of the information to be exported
    Return: none; outputs to a file named annotated_mutations.txt
    """
    print('Writing to file...')

    outputfile = open(name, 'w')

    #Output:
    #Write header
    outputfile.write('HGVS' + '\t'
                     + 'RSID' + '\t'
                     + 'VarType' + '\t'
                     + 'gnomADG' + '\t'
                     + 'gnomADE' + '\t'
                     + 'ClinVar' + '\t'
                     + 'Most Severe Consequence' + '\t'
                     + 'Genes Affected' + '\t'
                     + 'Relevant Transcripts' + '\n')

    #Write data
    i = 0
    #Write until but not including last element
    while i < (len(listanno) - 1):

        genelist = []
        relevanttranscripts = []

        #For genelist
        if listanno[i]['genelist'][0] is not None:
            for gene in listanno[i]['genelist']:
                 genelist.append(gene['gene_symbol'])
        #For relevanttranscripts
        if listanno[i]['relevanttranscripts'] is not None:
            for transcript in listanno[i]['relevanttranscripts']:
                relevanttranscripts.append(transcript['transcript_id'])

        outputfile.write(str(listanno[i]['_id']) + '\t'
                         + str(listanno[i]['rsid']) + '\t'
                         + str(listanno[i]['vartype']) + '\t'
                         + str(listanno[i]['gnomADG']) + '\t'
                         + str(listanno[i]['gnomADE']) + '\t'
                         + str(listanno[i]['ClinVar']) + '\t'
                         + str(listanno[i]['MScon']) + '\t'
                         + str(genelist) + '\t'
                         + str(relevanttranscripts) + '\t'
                         )

        if 'brain_expression' in listanno[i]:
            if listanno[i]['brain_expression'] is not None:
                outputfile.write(str(listanno[i]['brain_expression']))

        outputfile.write('\n')

        #Progress Indicator
        if i%1000 == 0:
            print (str(i) + ' out of ' + str(len(listanno)) + ' written...')

        i = i + 1
    #Write last element
    last_genelist = []
    last_relevanttranscripts = []

    if listanno[i]['genelist'][0] is not None:
        for gene in listanno[i]['genelist']:
             last_genelist.append(gene['gene_symbol'])
    if listanno[i]['relevanttranscripts'] is not None:
        for transcript in listanno[i]['relevanttranscripts']:
            last_relevanttranscripts.append(transcript['transcript_id'])

    outputfile.write(str(listanno[i]['_id']) + '\t'
                     + str(listanno[i]['rsid']) + '\t'
                     + str(listanno[i]['vartype']) + '\t'
                     + str(listanno[i]['gnomADG']) + '\t'
                     + str(listanno[i]['gnomADE']) + '\t'
                     + str(listanno[i]['ClinVar']) + '\t'
                     + str(listanno[i]['MScon']) + '\t'
                     + str(last_genelist) + '\t'
                     + str(last_relevanttranscripts) + '\t'
                     )

    if 'brain_expression' in listanno[i]:
        if listanno[i]['brain_expression'] is not None:
            outputfile.write(str(listanno[i]['brain_expression']))

    print('Annotations written to file')

    outputfile.close()
def outputsequence(sequence, id, type, path):
    sequencefile = open(path.as_posix()
                    + '/'
                    + id
                    + '_'
                    + type
                    + '.txt',
                    'w',)
    sequencefile.write(str(sequence))
    sequencefile.close()
def exportexcel(listcandidate):
    #Time
    now = datetime.datetime.now()
    time = now.strftime('%y%m%d-%H%M%S')
    #Job ID
    file_ID = input('Enter a job ID: ')

    #Create workbook
    wb = Workbook()
    filename = 'candidate_' + file_ID + '_' + time + '.xlsx'
    ws = wb.active
    ws.title = 'Results'

    #Write the header


    #Save
    wb.save(filename = filename)

#Parsing functions
def dumpvartype(idHGVS):
    """
    dumpvartype - determine the type of variant based on HGVS id
    Parameters: idHGVS - the HGVS id of variant being considered
    Return: vartype
    """
    vartype = ''
    if '>' in idHGVS:
        vartype = 'snv'
    elif 'del' in idHGVS:
        vartype = 'del'
    elif 'ins' in idHGVS:
        vartype = 'ins'
    return vartype
def converthgvs(idHGVS):
    """
    converthgvs - convert a HGVS ID into a format that can be used in a filename
    or directory
    Parameters: idHGVS - the HGVS ID as a string
    Return: a string with the converted HGVS ID
    """
    def findchr(idHGVS):
        if idHGVS[3:5].isnumeric() is True:
            return idHGVS[3:5]
        else:
            return idHGVS[3]
    def findpos(idHGVS, marker):
        pos_start = idHGVS.index(marker) + 1
        pos_end = pos_start + 1
        pos_complete = False

        while pos_complete is False:
            if idHGVS[pos_start:pos_end].isnumeric() is False:
                pos_chr = idHGVS[pos_start:(pos_end-1)]
                pos_complete = True
            else:
                pos_end = pos_end + 1
        return pos_chr
    #For substitution
    if '>' in idHGVS:
        #Identify chromosome number
        num_chr = findchr(idHGVS)
        #Identify position
        pos_chr = findpos(idHGVS, '.')
        #Identify allele 1 and 2
        allele1 = idHGVS[idHGVS.index('>') - 1]
        allele2 = idHGVS[idHGVS.index('>') + 1]
        return (num_chr + '-' + pos_chr + '-' + allele1 + '-' + allele2)

    #For deletion and insertion
    elif 'delins' in idHGVS:
        #Identify chromosome numbers
        num_chr = findchr(idHGVS)
        #Identify position
        pos_chr = findpos(idHGVS, '.')
        #Identify inserted sequence
        seq_delins = idHGVS[idHGVS.index('s') + 1:]
        #Debug
        return(num_chr + '-' + pos_chr + '-delins-' + seq_delins)

    #For deletion
    elif 'del' in idHGVS:
        #Identify chromosome number
        num_chr = findchr(idHGVS)
        #Identify position 1
        if '_' in idHGVS:
            pos1_chr = idHGVS[(idHGVS.index('.') + 1):idHGVS.index('_')]
            #Identify position 2
            pos2_chr = findpos(idHGVS, '_')
            return (num_chr + '-' + str(pos1_chr) + '_' + str(pos2_chr) + '-del')
        else:
            pos_chr = findpos(idHGVS, '.')
            return (num_chr + '-' + str(pos_chr) + '-del')

    #For insertion
    elif 'ins' in idHGVS:
        #Identify chromosome number
        num_chr = findchr(idHGVS)
        #Identify positions
        pos1_chr = idHGVS[(idHGVS.index('.') + 1):idHGVS.index('_')]
        pos2_chr = findpos(idHGVS, '_')
        #Identify inserted sequences
        seq_ins = idHGVS[idHGVS.index('s') + 1:]
        return (num_chr + '-' + str(pos1_chr) + '_' + str(pos2_chr) + '-ins' + seq_ins)
#Filtering functions
def filteraffected(listaffected, listcontrol):
    """
    filteraffected - filters from lists of HGVS IDs of affected and control
                     individuals
    Parameters:
    listaffected - list - contains the list of HGVS ids of affected individuals
    listcontrol - list - contains the list of HGVS ids of unaffected individuals
    Returns: list containing candidate mutations only
    """
    i = 0
    listfiltered = []

    print('Filtering starting...')

    while i < len(listaffected[0]):
        candidate = True

        #Check if HGVS is in all listAffected. If not, set candidate to false
        for a in listaffected:
            if (listaffected[0][i] in a) == False:
                candidate = False

        #Check if HGVS is in listControl. If yes, set candidate to False
        for c in listcontrol:
            if listaffected[0][i] in c:
                candidate = False

        #Write to file
        if candidate == True:
            listfiltered.append(listaffected[0][i].strip('\n'))

        if i%1000 == 0:
            print(str(i) + '/' + str(len(listaffected[0])) + ' completed...')

        i = i + 1

    print('Filtering completed.')
    return listfiltered
def filternodata(listanno):
    list_nodata = []
    list_candidate = []

    for anno in listanno:
        nodata = checknodata(anno)
        if nodata is True:
            list_nodata.append(anno)
        elif nodata is False:
            list_candidate.append(anno)

    return list_nodata, list_candidate
def filterfreq(listanno):
    list_overfreq = []
    list_candidate = []

    for anno in listanno:
        overfreq = checkfreq(anno)
        if overfreq is True:
            list_overfreq.append(anno)
        elif overfreq is False:
            list_candidate.append(anno)

    return list_overfreq, list_candidate
def filtercons(listanno):
    list_irrelevant = []
    list_candidate = []

    for anno in listanno:
        irrelevant = checkcons(anno)
        if irrelevant is True:
            list_irrelevant.append(anno)
        elif irrelevant is False:
            list_candidate.append(anno)

    return list_irrelevant, list_candidate
def filterexpression(listanno):
    list_hpa = []
    list_notbrainexpressed = []
    list_candidate = []

    #Get annotation from Human Protein Atlas
    for anno in listanno:
        if anno['genes'] is not None:
            for gene in anno['genes']:
                anno_hpa = hpa.annotate(gene['gene_id'])
                list_hpa.append(anno_hpa)
                if anno_hpa is not None:
                    gene['RNAbrd'] = anno_hpa['RNA brain regional distribution']
                else:
                    gene['RNAbrd'] = None

    #Filtering
    for anno in listanno:
        notbrainexpressflag = True
        if anno['genes'] is not None:
            for gene in anno['genes']:
                if gene['RNAbrd'] != 'Not detected':
                    notbrainexpressflag = False

        if notbrainexpressflag is True:
            list_notbrainexpressed.append(anno)
        else:
            list_candidate.append(anno)

    #Export and return
    exportanno(list_hpa, 'hpa_annotations.txt')
    return list_notbrainexpressed, list_candidate
def checknodata(anno):
    """
    checknodata - checks if variant has any data
    Parameters: data - variant and its associated annotations
    Returns: True or False
    """
    nodataflag = False
    if anno['MScon'] == 'N/A':
        nodataflag = True

    return nodataflag
def checkfreq(anno):
    """
    checkfreq - checks if variant above the frequency cut off
    Parameters: data - variant and its associated annotations
    Returns: True or False
    """

    overfreqpcflag = False
    if (anno['gnomADG'] != 'N/A') and (anno['gnomADG'] is not None):
        if anno['gnomADG'] >= 0.01:
            overfreqpcflag = True
    if (anno['gnomADE'] != 'N/A') and (anno['gnomADE'] != None):
        if anno['gnomADE'] >= 0.01:
            overfreqpcflag = True

    return overfreqpcflag
def checkcons(anno):
    """
    checkcons - checks if the variant is relevant
    Parameters: data - variant and its associated annotations
    Returns: True or False
    """
    nonrelevantflag = False

    if (anno['MScon'] == 'intron_variant'
        or anno['MScon'] == 'non_coding_transcript_exon_variant'
        or anno['MScon'] == 'upstream_gene_variant'
        or anno['MScon'] == 'downstream_gene_variant'
        or anno['MScon'] == 'synonymous_variant'
        or anno['MScon'] == '5_prime_UTR_variant'
        or anno['MScon'] == '3_prime_UTR_variant'
        or anno['MScon'] == 'intergenic_variant'):
        nonrelevantflag = True

    return nonrelevantflag
#Processing functions
def combineanno(listmvi, listvep, listHGVS):
    """
    combineanno - combine annotations from myvariant.info and VEP
    Parameters:
    listmvi - list containing dictionaries of myvariant.info annotations
    listvep - list containing dictionaries of VEP annotations
    listHGVS - list containing HGVS IDs
    Return: return combined annotation list
    """
    #listanno - holds list of combined annotations
    listanno = []
    i = 0

    #Do for every variant in listHGVS
    while i < len(listHGVS):
        #True if there's data
        anno_mvi_flag = False
        anno_vep_flag = False

        #Check if there's data
        if listmvi[i] is not None:
            anno_mvi_flag = True
        if listvep[i] is not None:
            anno_vep_flag = True

        #Progress tracker
        if i%100 == 0:
            print(str(i) + ' out of ' + str(len(listHGVS)) + ' completed..')

        #Combine annotations
        listanno.append(dict({'_id':listHGVS[i],
                        'rsid':mvi.dumpRSID(listmvi[i]),
                        'vartype':dumpvartype(listHGVS[i]),
                        'gnomADG':mvi.dumpgnomADG(listmvi[i]),
                        'gnomADE':mvi.dumpgnomADE(listmvi[i]),
                        'ClinVar':mvi.dumpCV(listmvi[i]),
                        'MScon':vep.dumpconsequence(listvep[i]),
                        'genes':vep.dumpensembltranscripts(listvep[i])
                        }))
        i = i + 1

    return listanno
def getsequences(listanno, basepath):
    def substitution(transcript, pos, allele, strand):
        if strand == -1:
            allele = str(Seq(allele, IUPAC.unambiguous_dna).complement())
        return transcript[:pos-1] + allele + transcript[pos:]

    def insertion(transcript, start, end, insert, strand):
        if strand == -1:
            insert = str(Seq(insert, IUPAC.unambiguous_dna).complement())
        return transcript[:start] + insert + transcript[end:]

    def deletion(transcript, start, end):
        return transcript[:start-1] + transcript[end:]

    def translation (sequence):
        if 'N' not in sequence:
            dna = Seq(sequence, IUPAC.unambiguous_dna)
            mrna = dna.transcribe()
            protein = mrna.translate()
            return str(protein)
        else:
            print('Ambiguous nucleotide present, no protein')

    #Do for every variant
    for anno in listanno:
        #Holds stats
        list_stats = []
        #Make a folder for variant
        print ('Working on ' + anno['_id'] + '...')
        variantpath = basepath.joinpath(converthgvs(anno['_id']))
        print(variantpath)
        variantpath.mkdir(exist_ok=True)
        #Do for every gene
        for gene in anno['genes']:
            #Stats to track
            list_transcripts = []
            list_notranscripts = []
            #Make folder for gene
            genepath = variantpath.joinpath(gene['gene_id'])
            genepath.mkdir(exist_ok=True)
            #Do for every transcript
            for transcript in gene['transcripts']:
                list_transcripts.append(transcript['transcript_id'])
                print('Pulling transcripts for ' + transcript['transcript_id'])
                cdna = vep.ensemblsequence(transcript['transcript_id'], 'cdna')
                cds = vep.ensemblsequence(transcript['transcript_id'], 'cds')
                #Write cdna and cds
                outputsequence(cdna, transcript['transcript_id'], 'cdna', genepath)
                outputsequence(cds, transcript['transcript_id'], 'cds', genepath)
                #Get rid of splice_region_variant, TF_binding_site, regulatory_region_variant
                if cds is not None:
                    get_protein = True
                else:
                    get_protein = False
                    list_notranscripts.append(transcript['transcript_id'])
                for term in transcript['consequence_terms']:
                    if (term == 'splice_region_variant'
                        or term == 'TF_binding_site'
                        or term == 'regulatory_region_variant'):
                        get_protein = False
                #Get protein sequences
                if (get_protein is True) and (anno['vartype'] == 'snv'):
                    vcds = substitution(cds,
                                        transcript['cds_start'],
                                        transcript['variant_allele'],
                                        transcript['strand'],
                                        )
                    vprotein = translation(vcds)
                elif (get_protein is True) and (anno['vartype'] == 'del'):
                    vcds = deletion(cds,
                                    transcript['cds_start'],
                                    transcript['cds_end'],
                                    )
                    vprotein = translation(vcds)
                elif (get_protein is True) and (anno['vartype'] == 'ins' or anno['vartype' == 'delins']):
                    vcds = insertion(cds,
                                    transcript['cds_start'],
                                    transcript['cds_end'],
                                    transcript['variant_allele'],
                                    transcript['strand']
                                    )
                    vprotein = translation(vcds)
                #Write proteins
                if get_protein is True:
                    protein = translation(cds)
                    outputsequence(protein, transcript['transcript_id'], 'protein', genepath)
                    outputsequence(vprotein, transcript['transcript_id'], 'vprotein', genepath)
            #Add to the stats document
            list_stats.append(dict({'gene':gene['gene_id'],
                                    'transcripts':list_transcripts,
                                    'no_transcripts':list_notranscripts,
                                    }))
        #Write the summary file for the variant
        summary_file = open(variantpath.as_posix()
                            + '/'
                            + converthgvs(anno['_id'])
                            + '_summary.txt'
                            ,
                            'w',
                            )
        for item in list_stats:
            summary_file.write(item['gene']
                                + '\n'
                                + 'All transcripts:'
                                + str(item['transcripts'])
                                + '\n'
                                + 'No sequence:'
                                + str(item['no_transcripts'])
                                + '\n'
                                )

        summary_file.close()

##### Program Start ##################################################################
#GLOBAL VARIABLES
EXIT_PROGRAM = False
#LIST_ANNO - holds all annotations
LIST_ANNO = []
#LIST_VEP - holds VEP annotations
LIST_VEP = []
#LIST_MVI - holds myvariant.info annotations
LIST_MVI = []
#LIST_CANDIDATE = holds the candidate annotations
LIST_CANDIDATE = []

BASEPATH = pathlib.Path.cwd()

while EXIT_PROGRAM == False:
    #Display menu
    print('1) Import HGVS from .vcf files')
    print('2) Generate filtered HGVS file')
    print('3) Annotation')
    print('4) Import data')
    print('5) Export data')
    print('6) Process data')
    print('7) Get Ensembl CDS sequences')
    print('8) Exit (modules still loaded)')

    #Accept input from user
    option = int(input('Select option: '))
    #Input check
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

        filename = input('Please enter the .vcf file to extract HGVS from: ')
        vcftoHGVS(filename)

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
        listfiltered= filteraffected(listaffected, listcontrol)
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
        end = time.time()
        print('Total run time: ' + str(end - start) + '\n')

    #3) Annotation
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
                filename = input('Please enter the name of the file listing HGVS IDs to be annotated: ')
                listHGVS = importHGVS(filename)
                #Perform MVI annotation
                start = time.time()
                LIST_MVI = mvi.annotate(listHGVS)
                end = time.time()
                print('Total run time: ' + str(end - start) + '\n')

            elif optionannotation == 2:
                #Open HGVS ID file that will be annotated
                filename = input('Please enter the name of the file listing HGVS IDs to be annotated: ')
                listHGVS = importHGVS(filename)
                #Perform VEP annotation
                start = time.time()
                LIST_VEP = vep.annotate(listHGVS)
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
                #Open the HGVS ID file
                filename = input('Please enter the name of the file being imported : ')
                LIST_MVI = importanno(filename)
                print('Data imported.' + '\n')
            elif optionimport == 2:
                filename = input('Please enter the name of the file being imported : ')
                LIST_VEP = importanno(filename)
                print('Data imported.' + '\n')
            elif optionimport == 3:
                filename = input('Please enter the name of the file being imported : ')
                LIST_ANNO = importanno(filename)
                print('Data imported.' + '\n')
            elif optionimport == 4:
                exitimport = True
                print('')

    elif option == 5:
        exitexport = False
        while exitexport == False:
            print('1) Export myvariant.info annotation data')
            print('2) Export VEP annotation data')
            print('3) Export combined data')
            print('4) Return to main menu')

            optionexport = int(input('Select option: '))

            while (optionexport != 1
                   and optionexport != 2
                   and optionexport != 3
                   and optionexport != 4):
                optionexport = int(input('Invalid selection; please select one of the options: '))

            if optionexport == 1:
                exportanno(LIST_MVI, 'myvariantannotation.txt')
            elif optionexport == 2:
                exportanno(LIST_VEP, 'vepannotation.txt')
            elif optionexport == 3:
                filename = input('Please enter the name of the file being imported : ')
                listHGVS = importHGVS(filename)
                LIST_ANNO = combineanno(LIST_MVI, LIST_VEP, listHGVS)
                exportanno(LIST_ANNO, 'annotated_mutations.txt')
            elif optionexport == 4:
                exitexport = True
                print('')

    elif option == 6:
        #Lists
        list_nodata = []
        list_irrelevant = []
        list_highfreq = []
        list_notexpressedbrain = []
        list_candidate = []
        list_intermediate = []

        #Job Identification
        now = datetime.datetime.now()
        time = now.strftime('%y%m%d-%H%M%S')

        #File-handling variables
        basepath = pathlib.Path.cwd()
        file_created = False
        while file_created is False:
            file_ID = input('Enter a job ID: ')
            newpath = basepath.joinpath(file_ID+'_'+time)
            try:
                newpath.mkdir()
            except FileExistsError as fee_error:
                print(fee_error)
                print('Please try another ID: ')
            else:
                file_created = True

        #Get combine annotations
        #Open HGVS ID file that will be annotated
        filename = input('Please enter the name of the file listing HGVS IDs to be annotated: ')
        listHGVS = importHGVS(filename)
        LIST_ANNO = combineanno(LIST_MVI, LIST_VEP, listHGVS)

        #Filenames
        filename_nodata = 'nodata_' + str(file_ID) + '_' + time + '.txt'
        filename_irrelevant = 'nonrelevant_' + str(file_ID) + '_' + time + '.txt'
        filename_highfreq = 'overfreqpc_' + str(file_ID) + '_' + time + '.txt'
        filename_notexpressedbrain = 'notbrainexpress_' + str(file_ID) + '_' + time + '.txt'
        filename_candidate = 'candidatevariants_' + str(file_ID) + '_' + time + '.txt'
        filename_summary = 'summary_' + str(file_ID) + '_' + time + '.txt'

        #Conduct filtering
        #No data
        print('Filtering out no data...')
        list_nodata, list_intermediate = filternodata(LIST_ANNO)
        print('No data variants filtered out.')
        #Consequence filter
        print('Filtering out irrelevant variants...')
        list_irrelevant, list_intermediate = filtercons(list_intermediate)
        print('Irrelevant variants filtered out.')
        #Frequency filter
        print('Filtering out high frequency variants...')
        list_highfreq, list_intermediate = filterfreq(list_intermediate)
        print('High frequency variants filtered out.')
        #Expression filter
        print('Filtering out variants not expressed in the brain...')
        list_notexpressedbrain, list_candidate = filterexpression(list_intermediate)
        print('Variants not expressed in brain filtered out.')

        LIST_CANDIDATE = list_candidate

        #Export all of the files
        exportanno(list_nodata, newpath.as_posix() + '/' + filename_nodata)
        exportanno(list_irrelevant, newpath.as_posix() + '/' + filename_irrelevant)
        exportanno(list_highfreq, newpath.as_posix() + '/' + filename_highfreq)
        exportanno(list_notexpressedbrain, newpath.as_posix() + '/' + filename_notexpressedbrain)
        exportanno(list_candidate, newpath.as_posix() + '/' + filename_candidate)

        #Export candidate by snp type
        list_snv = []
        list_del = []
        list_ins = []
        for candidate in list_candidate:
            if candidate['vartype'] == 'snv':
                print('snp')
                list_snv.append(candidate)
            elif candidate['vartype'] == 'del':
                print ('del')
                list_del.append(candidate)
            elif candidate['vartype'] == 'ins' or candidate['vartype'] == 'delins':
                print ('delins')
                list_ins.append(candidate)

        exportanno(list_snv,
                        newpath.as_posix()
                        + '/'
                        + 'candidatesnv_'
                        + str(file_ID)
                        + '_'
                        + time
                        + '.txt'
                        )
        exportanno(list_del,
                        newpath.as_posix()
                        + '/'
                        + 'candidatedel_'
                        + str(file_ID)
                        + '_'
                        + time
                        + '.txt'
                        )
        exportanno(list_ins,
                        newpath.as_posix()
                        + '/'
                        + 'candidateins_'
                        + str(file_ID)
                        + '_'
                        + time
                        + '.txt'
                        )
        #Summary
        summaryfile = open(newpath.as_posix() + '/' + filename_summary, 'w')

        summaryfile.write('Summary report of analysis'
                        + '\n'
                        + 'ID: '
                        + file_ID
                        + '\n'
                        + 'Date performed: '
                        + time
                        + '\n'
                        + 'Total number of samples: '
                        + str(len(LIST_ANNO))
                        + '\n'
                        + '# No Data: '
                        + str(len(list_nodata))
                        + '\n'
                        + '# Non-relevant: '
                        + str(len(list_irrelevant))
                        + '\n'
                        + '# Greater than 0.1% AF: '
                        + str(len(list_highfreq))
                        + '\n'
                        + '# Not expressed in brain: '
                        + str(len(list_notexpressedbrain))
                        + '\n'
                        + '# Candidates: '
                        + str(len(list_candidate))
                        + '\n'
                        + '# Candidate SNV:'
                        + str(len(list_snv))
                        + '\n'
                        + '# Candidate deletions:'
                        + str(len(list_del))
                        + '\n'
                        + '# Candidate insertions:'
                        + str(len(list_ins))
                        )
        summaryfile.close()
    elif option == 7:
        #Import annotations if not already in program:
        if len(LIST_CANDIDATE) == 0:
            filename = input('Please enter the name of the file being imported : ')
            LIST_CANDIDATE = importanno(filename)
        getsequences(LIST_CANDIDATE, BASEPATH)
    elif option == 8:
        check = input('Are you sure you want to exit? (y/n)')
        if check == 'y':
            EXIT_PROGRAM = True
        else:
            print('Returning to menu' + '\n')

print ("Terminating program...")
