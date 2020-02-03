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

    return listHGVS

def importanno(filename):
    """
    importanno - from a file containing json.dumps annotation data, import it as
    a list of dictionaries
    Parameters: none
    Return: a list of dictionaries of annotation data
    """

    inputfile_open = False
    listdata = []

    #Open the file
    while inputfile_open == False:
        try:
            inputfile = open(filename, 'r')
            inputfileopen = True
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

    outputFile.close()

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
                     + 'MSCon' + '\t'
                     + 'GeneList' + '\n')

    #Write data
    for i in listanno:

        outputfile.write(i['_id'] + '\t'
                         + i['rsid'] + '\t'
                         + i['vartype'] + '\t'
                         + str(i['gnomADG']) + '\t'
                         + str(i['gnomADE']) + '\t'
                         + str(i['ClinVar']) + '\t'
                         + i['MScon'] + '\t'
                         + json.dumps(i['genelist'])
                         )

        if listanno.index(i) != (len(listanno)-1):
            outputfile.write('\n')

        #Progress Indicator
        if listanno.index(i)%1000 == 0:
            print (str(listanno.index(i)) + ' out of ' + str(len(listanno)) + ' written...')

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
                                            'gnomad_genome.af', 'gnomad_exome.af']))

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

##### Parsing Functions ########################################################
def dumpCV(anno_mvi):
    """
    dumpCV - pull ClinVar data from the JSON data structure
    Parameters: data - dictionary of the JSON data on the variant
    Return: list containing number of pathogenic, uncertain significance, and benign reports
    """
    #0 - pathogenic, 1 - uncertain significance, 2 - benign
    num_report = [0, 0, 0]
    chk_flag = False

    #Check if 'clinvar' exists
    if 'clinvar' in anno_mvi:
        #Check if 'rcv' exists
        if 'rcv' in anno_mvi['clinvar']:
            num_rcv = len(anno_mvi['clinvar']['rcv'])
            #Number of RCV entries affects parsing of the json data
            if num_rcv > 1:
                for rcv in anno_mvi['clinvar']['rcv']:
                    if (rcv['clinical_significance'] == 'Pathogenic'
                        or rcv['clinical_significance'] == 'Likely pathogenic'
                        or rcv['clinical_significance'] == 'Pathogenic/Likely pathogenic'):
                        numReport[0] = numReport[0] + 1
                    elif rcv['clinical_significance'] == 'Uncertain significance':
                        numReport[1] = numReport[1] + 1
                    elif (rcv['clinical_significance'] == 'Benign'
                          or rcv['clinical_significance'] == 'Likely benign'
                          or rcv['clinical_significance'] == 'Benign/likely benign'):
                        numReport[2] = numReport[2] + 1
                    else:
                        chk_flag = True
            if num_rcv <= 1:
                if (anno_mvi['clinvar']['rcv']['clinical_significance'] == 'Pathogenic'
                    or anno_mvi['clinvar']['rcv']['clinical_significance'] == 'Likely pathogenic'
                    or anno_mvi['clinvar']['rcv']['clinical_significance'] =='Pathogenic/Likely pathogenic'):
                    numReport[0] = 1
                elif anno_mvi['clinvar']['rcv']['clinical_significance'] == 'Uncertain significance':
                    numReport[1] = 1
                elif (anno_mvi['clinvar']['rcv']['clinical_significance'] == 'Benign'
                      or anno_mvi['clinvar']['rcv']['clinical_significance'] == 'Likely benign'
                      or anno_mvi['clinvar']['rcv']['clinical_significance'] =='Benign/Likely benign'):
                    numReport[2] = 1
                else:
                    chk_flag = True

            return num_report

    else:
        return None

def dumpRSID(anno_mvi):
    """
    dumpRSID - pull the dbSNP RSID from the JSON data structure
    Parameters: data - dictionary/list of the JSON data on the variant
    Return: the RSID
    """
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
    Parameters: data - dictionary/list of the JSON data on the variant
    Return: genome allele frequency
    """
    if 'gnomad_genome' in anno_mvi:
        if 'af' in anno_mvi['gnomad_genome']:
            if 'af' in anno_mvi['gnomad_genome']['af']:
                return anno_mvi['gnomad_genome']['af']['af']
    else:
        return None

def dumpgnomADE(anno_mvi):
    """
    dumpgnomADE - pull the gnomAD allele frequency (excl. genome sequences) from the JSON data structure
    Parameters: data - dictionary/list of the JSON data on the variant
    Return: exome allele frequency
    """
    if 'gnomad_exome' in anno_mvi:
        if 'af' in anno_mvi['gnomad_exome']:
            if 'af' in anno_mvi['gnomad_exome']['af']:
                return anno_mvi['gnomad_exome']['af']['af']
    else:
        return None

def dumpconsequence(anno_vep):
    """
    dumpconsequence - pull the most severe consequence from the JSON data structure
    Parameters: data - dictionary/list of the JSON data on the variant
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
    Parameters: data - dictionary/list of the JSON data on the variant
    Return: list of gene IDs that the variant affects
    """

    geneids = []
    i = 0

    if anno_vep is not None:
        if 'transcript_consequences' in anno_vep:
            while i < len(anno_vep['transcript_consequences']):
                #Check if the gene_id is in the current list
                names = {'gene_id': anno_vep['transcript_consequences'][i]['gene_id'],
                        'gene_symbol': anno_vep['transcript_consequences'][i]['gene_symbol']
                        }
                if not(names in geneids):
                    #geneids.append(data['transcript_consequences'][i]['gene_id'])
                    geneids.append(names)

                i = i + 1
        else:
            geneids.append(None)
    else:
        geneids.append(None)

    return geneids

def combineanno(listmvi, listvep):
    """

    combineanno - combines annotations from MVI and VEP together into one list of dictionaries
    Parameters:
    listmvi - list - contains a list of dictionaries of MVI information
    listvep - list - contains a list of dictionaries of VEP information
    Returns: a list containing both sets of information

    """
    al = []
    lc = importHGVS()
    i = 0

    while i < len(lc):
        if i%100 == 0:
            print(str(i) + ' out of ' + str(len(lc)) + ' completed..')
        if str(type(listmvi[i])) != "<class 'NoneType'>":
            al.append(dict({'_id':listmvi[i]['_id'],
                            'rsid':dumpRSID(listmvi[i]),
                            'vartype':dumpvartype(listmvi[i]),
                            'gnomADG':dumpgnomADG(listmvi[i]),
                            'gnomADE':dumpgnomADE(listmvi[i]),
                            'ClinVar':dumpCV(listmvi[i]),
                            'MScon':dumpconsequence(listvep[i]),
                            'genelist':dumpensemblgeneid(listvep[i])
                            }))

        if str(type(listmvi[i])) == "<class 'NoneType'>":
            al.append(dict({'_id':lc[i],
                            'rsid':'N/A',
                            'vartype':'N/A',
                            'gnomADG':'N/A',
                            'gnomADE':'N/A',
                            'ClinVar':'N/A',
                            'MScon':dumpconsequence(listvep[i]),
                            'genelist':dumpensemblgeneid(listvep[i])
                            }))

        i = i + 1

    return al

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

def filternodata(anno):
    """
    filternodata - checks if variant has any data
    Parameters: data - variant and its associated annotations
    Returns: True or False
    """
    nodataflag = False
    if anno['MScon'] == 'N/A':
        nodataflag = True

    return nodataflag

def filterfreq(anno):
    """
    filterfreq - checks if variant above the frequency cut off
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

def filtercons(anno):
    """
    filtercons - checks if the variant is relevant
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

##### Move to VariantAnnotation.py #############################################
def filtervariant(listanno, name = ""):
    """
    filtervariant - filters variants based on selected criteria
    Parameters: listanno - list of the annotated variants
    Returns: returns the filtered list of annotated variants
    """
    start = time.time()

    nodata = []
    nonrelevant = []
    overfreqpc = []
    filterstep3 = []

    notbrainexpress = []
    expressiondata = []

    candidate = []
    now = datetime.datetime.now()

    if name == "":
        name = now.strftime('%y%m%d-%H%M%S')

    for data in listanno:
        #Filter flags
        nodataflag = filternodata(data)
        nonrelevantflag = filtercons(data)
        overfreqpcflag = filterfreq(data)

        #Append data depending on which flags were triggered. If none were triggered
        #then append to candidate list
        if nodataflag == True:
            nodata.append(data)
        elif nonrelevantflag == True:
            nonrelevant.append(data)
        elif overfreqpcflag == True:
            overfreqpc.append(data)
        else:
            filterstep3.append(data)
        if listanno.index(data)%1000 == 0:
            print (str(listanno.index(data)) + ' out of ' + str(len(listanno)) + ' filtered...')

    #HPA annotation
    for data in filterstep3:

        print(str(filterstep3.index(data)))
        print(data['_id'])

        variantexpression = []
        numnotdetected = 0
        notbrainexpressflag = False

        if data['genelist'][0] != 'N/A':
            x = annothpa(data)
            for gene in x:
                if gene['RNA brain regional distribution'] == 'Not detected':
                    #notbrainexpressflag = True
                    numnotdetected = numnotdetected + 1
                variantexpression.append(dict({
                                            'gene':gene['Gene'],
                                            'RNAbrd':gene['RNA brain regional distribution']
                }))
            if numnotdetected == len(data['genelist']):
                notbrainexpressflag = True
        if data['genelist'][0] == 'N/A':
            x = None
            variantexpression.append(None)
        expressiondata.append(x)

        data['brain_expression'] = variantexpression

        if notbrainexpressflag == True:
            notbrainexpress.append(data)
        else:
            candidate.append(data)

    exportanno(expressiondata, 'expressiondata.txt')

    #Write to files
    writeanno(nodata, 'nodata_'+ name + '.txt')
    writeanno(nonrelevant, 'nonrelevant_' + name + '.txt')
    writeanno(overfreqpc, 'overfreqpc_' + name + '.txt')

    #Output not brainexpress:
    outputfilenbe = open('notbrainexpress_' + name + '.txt', 'w')
    #Write header
    outputfilenbe.write('HGVS' + '\t'
                     + 'RSID' + '\t'
                     + 'VarType' + '\t'
                     + 'gnomADG' + '\t'
                     + 'gnomADE' + '\t'
                     + 'ClinVar' + '\t'
                     + 'MSCon' + '\t'
                     + 'GeneList' + '\t'
                     + 'BrainExpression' + '\n')
    #Write data
    for data in notbrainexpress:
        outputfilenbe.write(data['_id'] + '\t'
                         + data['rsid'] + '\t'
                         + data['vartype'] + '\t'
                         + str(data['gnomADG']) + '\t'
                         + str(data['gnomADE']) + '\t'
                         + str(data['ClinVar']) + '\t'
                         + data['MScon'] + '\t'
                         + json.dumps(data['genelist']) + '\t'
                         + json.dumps(data['brain_expression']) + '\n'
                         )
    outputfilenbe.close()

    #Output candidate:
    outputfilec = open('candidates_' + name + '.txt', 'w')
    #Write header
    outputfilec.write('HGVS' + '\t'
                     + 'RSID' + '\t'
                     + 'VarType' + '\t'
                     + 'gnomADG' + '\t'
                     + 'gnomADE' + '\t'
                     + 'ClinVar' + '\t'
                     + 'MSCon' + '\t'
                     + 'GeneList' + '\t'
                     + 'BrainExpression' + '\n')
    #Write data
    for data in candidate:
        outputfilec.write(data['_id'] + '\t'
                         + data['rsid'] + '\t'
                         + data['vartype'] + '\t'
                         + str(data['gnomADG']) + '\t'
                         + str(data['gnomADE']) + '\t'
                         + str(data['ClinVar']) + '\t'
                         + data['MScon'] + '\t'
                         + json.dumps(data['genelist']) + '\t'
                         + json.dumps(data['brain_expression']) + '\n'
                         )
    outputfilenbe.close()

    end = time.time()
    print('Processing time: ' + str(end - start))

    #Summary report
    identifier = input('Please enter a Job ID (description): ')
    outputsummary = open('summary_' + name + '.txt', 'w')
    outputsummary.write('Summary report of analysis'
                        + '\n'
                        + 'ID: '
                        + identifier
                        + '\n'
                        + 'Date performed: '
                        + name
                        + '\n'
                        + 'Total number of samples: '
                        + str(len(listanno))
                        + '\n'
                        + '# No Data: '
                        + str(len(nodata))
                        + '\n'
                        + '# Non-relevant: '
                        + str(len(nonrelevant))
                        + '\n'
                        + '# Greater than 0.1% AF: '
                        + str(len(overfreqpc))
                        + '\n'
                        + '# Not expressed in brain: '
                        + str(len(notbrainexpress))
                        + '\n'
                        + '# Candidates: '
                        + str(len(candidate))
                        + '\n'
                        + 'Processing time: '
                        + str (end - start)
                        + ' seconds'
                        )

    return candidate
