"""
va.py - wrapper for myvariant.info and Ensembl's VEP

"""
##### Function definitions #####################################################
import json
import time
import requests
import sys
import myvariant
import datetime

##### Importing Functions ######################################################
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

##### Exporting Functions ######################################################
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

##### Annotation functions #####################################################
def annotmvi(listHGVS):
    """
    annotmvi - accepts HGVS, returns data on mutation
    Parameters: listHGVS: list of HGVS IDs to retrieve annotations for
    Return: dictionary/list containing the json data from myvariant.info
    """
    listmvi = []
    mv = myvariant.MyVariantInfo()

    #For each HGVS ID in the list, retrieve the annotation
    for idHGVS in listHGVS:
        listmvi.append(mv.getvariant(idHGVS, fields = ['dbsnp.rsid',
                                            'dbsnp.alleles',
                                            'dbsnp.vartype',
                                            'dbsnp.gene',
                                            'clinvar',
                                            'gnomad_genome.af',
                                            'gnomad_exome.af',
                                            'dbnsfp.ensembl',
                                            'dbnsfp.uniprot',
                                            'dbnsfp.polyphen2',
                                            'dbnsfp.sift',
                                            'dbnsfp.provean',
                                            ]))

        #Progress tracker
        if listHGVS.index(idHGVS)%100 == 0:
            print (str(listHGVS.index(idHGVS))
                    + ' out of '
                    + str(len(listHGVS))
                    + ' written...')

    return listmvi

def annotvep(listHGVS):
    """
    annotvep - from a list of HGVS IDs, retrieve all VEP annotations
    Paramters:
    lc - list - list of HGVS ids
    Return: a list of dictionaries of VEP data
    """
    #Processing time tracker
    start = time.time()

    #indexes for list traversal
    i = 0
    u = 200

    #VEP REST API
    server = "http://grch37.rest.ensembl.org"
    ext = "/vep/human/hgvs"
    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}

    #listvep - will store all of the data
    listvep = []

    #Annotation process
    print('Initiating annotation...')
    while i < len(listHGVS):

        #Retrieve 200 HGVS IDs - maximum for the batch query to the VEP REST API
        if u != (len(listHGVS) - 1):
            rangeHGVS = listHGVS[i:u]
        if u == (len(listHGVS) - 1):
            rangeHGVS = listHGVS[i:]
        p_rangeHGVS = str(rangeHGVS).replace("'", '"')

        #Retrieval of data from VEP REST API
        r = requests.post(server+ext, headers=headers, data=('{ "hgvs_notations" : ' + p_rangeHGVS + ' }'))

        #Error handling
        if not r.ok:
            r.raise_for_status()
            sys.exit()

        #Store as list of dictionaries
        decoded = r.json()

        #Check for missing annotations
        if len(rangeHGVS) != len(decoded):
            print(str(len(rangeHGVS)
                    - len(decoded))
                    + ' annotations are missing')
            try:
                q = 0
                for idHGVS in rangeHGVS:
                    if rangeHGVS.index(idHGVS) < len(decoded):
                        if idHGVS != decoded[q]['id']:
                            decoded.insert(rangeHGVS.index(idHGVS), None)
                        q = q + 1

                    if rangeHGVS.index(hgvs) >= len(decoded):
                        decoded.append(None)
            except:
                print('Issue between '
                        + str(rangeHGVS[0])
                        + 'and '
                        + str(rangeHGVS[len(rangeHGVS)-1]))
                print(len(decoded))

        #Add these results to annot list
        for anno in decoded:
            listvep.append(anno)

        print (str(u) + ' out of ' + str(len(listHGVS)) + ' completed...')

        i = i + 200
        u = u + 200

        if u > len(listHGVS):
            u = len(listHGVS) - 1

    #Processing time tracker
    end = time.time()
    print('Processing time: ' + str(end - start))

    return listvep

def annothpa(data):
    """
    annothpa - pulls expression data from the Human Protein Atlas
    Parameters: data - program-generated annotations for a single variant
    Returns: returns expression data
    """
    #expression - list that holds the return data from
    expression = []
    i = 0

    server = "http://www.proteinatlas.org/"

    #For every gene in the gene list

    while i < len(data['genelist']):

        r = requests.get(server + data['genelist'][i]['gene_id'] + '.json')

        #Check if there was data returned
        if r.status_code == requests.codes.ok:
            expression.append(r.json())

        i = i + 1

    return expression

def ensemblsequence(transcriptid, seq):
    server = 'https://rest.ensembl.org'
    ext = '/sequence/id/' + transcriptid + '?'

    xheaders={'Content-Type':'application/json'}
    xparams = {'type':seq,'object_type':'transcript'}

    r = requests.get(server+ext, headers=xheaders, params=xparams)

    if r.status_code == 200:
        return r.json()['seq']
    else:
        return None

##### Parsing Functions ########################################################
def dumpCV(anno_mvi):
    """
    dumpCV - pull ClinVar data from the JSON data structure
    Parameters: anno_mvi - dictionary of the JSON data on the variant
    Return: list containing number of pathogenic, uncertain significance, and benign reports
    """
    #0 - pathogenic, 1 - uncertain significance, 2 - benign
    num_report = [0, 0, 0]
    chk_flag = False

    #Check if there's an annotation
    if anno_mvi is not None:
        #Check if 'clinvar' exists
        if 'clinvar' in anno_mvi:
            #Check if 'rcv' exists
            if 'rcv' in anno_mvi['clinvar']:
                if isinstance(anno_mvi['clinvar']['rcv'], list):
                    for rcv in anno_mvi['clinvar']['rcv']:
                        if (rcv['clinical_significance'] == 'Pathogenic'
                            or rcv['clinical_significance'] == 'Likely pathogenic'
                            or rcv['clinical_significance'] == 'Pathogenic/Likely pathogenic'):
                            num_report[0] = num_report[0] + 1
                        elif rcv['clinical_significance'] == 'Uncertain significance':
                            num_report[1] = num_report[1] + 1
                        elif (rcv['clinical_significance'] == 'Benign'
                              or rcv['clinical_significance'] == 'Likely benign'
                              or rcv['clinical_significance'] == 'Benign/likely benign'):
                            num_report[2] = num_report[2] + 1
                        else:
                            chk_flag = True
                else:
                    if (anno_mvi['clinvar']['rcv']['clinical_significance'] == 'Pathogenic'
                        or anno_mvi['clinvar']['rcv']['clinical_significance'] == 'Likely pathogenic'
                        or anno_mvi['clinvar']['rcv']['clinical_significance'] =='Pathogenic/Likely pathogenic'):
                        num_report[0] = 1
                    elif anno_mvi['clinvar']['rcv']['clinical_significance'] == 'Uncertain significance':
                        num_report[1] = 1
                    elif (anno_mvi['clinvar']['rcv']['clinical_significance'] == 'Benign'
                          or anno_mvi['clinvar']['rcv']['clinical_significance'] == 'Likely benign'
                          or anno_mvi['clinvar']['rcv']['clinical_significance'] =='Benign/Likely benign'):
                        num_report[2] = 1
                    else:
                        chk_flag = True

                return num_report

    else:
        return None

def dumpRSID(anno_mvi):
    """
    dumpRSID - pull the dbSNP RSID from the JSON data structure
    Parameters: anno_mvi - dictionary/list of the JSON data on the variant
    Return: the RSID
    """
    if anno_mvi is not None:
        if 'dbsnp' in anno_mvi:
            return anno_mvi['dbsnp']['rsid']
    else:
        return None

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

def dumpgnomADG(anno_mvi):
    """
    dumpgnomADG - pull the gnomAD allele frequency (inc. genome sequences) from the JSON data structure
    Parameters: anno_mvi - dictionary/list of the JSON data on the variant
    Return: genome allele frequency
    """
    if anno_mvi is not None:
        if 'gnomad_genome' in anno_mvi:
            if 'af' in anno_mvi['gnomad_genome']:
                if 'af' in anno_mvi['gnomad_genome']['af']:
                    return anno_mvi['gnomad_genome']['af']['af']
    else:
        return None

def dumpgnomADE(anno_mvi):
    """
    dumpgnomADE - pull the gnomAD allele frequency (excl. genome sequences) from the JSON data structure
    Parameters: anno_mvi - dictionary/list of the JSON data on the variant
    Return: exome allele frequency
    """
    if anno_mvi is not None:
        if 'gnomad_exome' in anno_mvi:
            if 'af' in anno_mvi['gnomad_exome']:
                if 'af' in anno_mvi['gnomad_exome']['af']:
                    return anno_mvi['gnomad_exome']['af']['af']
    else:
        return None

def dumpconsequence(anno_vep):
    """
    dumpconsequence - pull the most severe consequence from the JSON data structure
    Parameters: anno_vep - dictionary/list of the JSON data on the variant
    Return: most severe consequence
    """

    consequence = ""

    if anno_vep is not None:
        consequence = anno_vep['most_severe_consequence']
    else:
        consequence = None
    return consequence

def dumpensemblgeneid(anno_vep):
    """
    dumpensemblgeneid - pull the gene IDs from the JSON data structure
    Parameters: anno_vep - dictionary/list of the JSON data on the variant
    Return: list of gene IDs that the variant affects
    """

    geneids = []
    i = 0
    relevant_transcripts = dumprelevanttranscripts(anno_vep)

    if relevant_transcripts != []:
        for transcript in relevant_transcripts:
            #Check if the gene_id is in the current list
            names = {'gene_id': transcript['gene_id'],
                    'gene_symbol': transcript['gene_symbol']
                    }

            if not(names in geneids):
                geneids.append(names)
    else:
        geneids.append(None)

    return geneids

def dumprelevanttranscripts(anno_vep):
    """
    dumprelevanttranscripts - pull Ensembl transcript IDs
    Parameters: anno_vep - dictionary/list of the JSON data on the variant
    Return: list of relevant transcipts
    """
    #list_transcripts - holds all of the relevant transcripts
    list_transcripts = []
    nodataflag = True

    #Check if the annotation exists
    if anno_vep is not None:
        #Check if 'transcript_consequences' exists
        if 'transcript_consequences' in anno_vep:
            nodataflag = False

    #Do for each transcript
    if nodataflag is False:
        for transcript in anno_vep['transcript_consequences']:
            relevantflag = False
            #Do for each consequence term
            for term in transcript['consequence_terms']:
                if (term != 'intron_variant'
                    and term != 'non_coding_transcript_exon_variant'
                    and term != 'non_coding_transcript_variant'
                    and term != 'upstream_gene_variant'
                    and term != 'downstream_gene_variant'
                    and term != 'synonymous_variant'
                    and term != '5_prime_UTR_variant'
                    and term != '3_prime_UTR_variant'
                    and term != 'NMD_transcript_variant'):
                    relevantflag = True

            #Append to the list f relevant transcripts if it is relevant
            if relevantflag is True:
                list_transcripts.append(transcript)

    return list_transcripts

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
                        'rsid':dumpRSID(listmvi[i]),
                        'vartype':dumpvartype(listHGVS[i]),
                        'gnomADG':dumpgnomADG(listmvi[i]),
                        'gnomADE':dumpgnomADE(listmvi[i]),
                        'ClinVar':dumpCV(listmvi[i]),
                        'MScon':dumpconsequence(listvep[i]),
                        'genelist':dumpensemblgeneid(listvep[i]),
                        'relevanttranscripts':dumprelevanttranscripts(listvep[i])
                        }))
        i = i + 1

    return listanno

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

##### Processing Functions #####################################################
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
        if anno['gnomADG'] >= 0.001:
            overfreqpcflag = True
    if (anno['gnomADE'] != 'N/A') and (anno['gnomADE'] != None):
        if anno['gnomADE'] >= 0.001:
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
