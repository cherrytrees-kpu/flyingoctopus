"""
va.py - wrapper for myvariant.info and Ensembl's VEP

"""
##### Function definitions #######################################################################
import json
import time
import requests
import sys
import myvariant
import datetime

def outputHGVS(listHGVS, name = ""):
    """

    outputHGVS - outputs a file with the name 'HGVS_*name.txt' containing the list of HGVS ids
    Parameters:
    listHGVS - list - contains list ids
    Returns: nothing; outputs to file 'HGVS_*name.txt'

    """
    index = 0

    #Accept name input from user
    if name != "":
        name = input('Please enter an identifier for your file: ')
    outputFile = open('HGVS_' + name + '.txt', 'w')

    #Write to file
    while(len(listHGVS)) > index:
        outputFile.write(listHGVS[index])
        if index != (len(listHGVS)-1):
            outputFile.write('\n')
        index = index + 1

    outputFile.close()

def filteraffected(listAffected, listControl):
    """
    filteraffected - filters from lists of HGVS IDs of affected and control individuals
    Parameters:
    listAffected - list - contains the list of HGVS ids of affected individuals
    listControl - list - contains the list of HGVS ids of unaffected individuals
    Returns: list containing candidate mutations only
    """
    i = 0
    listCandidate = []

    print('Filtering starting...')

    while i < len(listAffected[0]):
        candidate = True

        #Check if HGVS is in all listAffected. If not, set candidate to false
        for a in listAffected:
            if (listAffected[0][i] in a) == False:
                candidate = False

        #Check if HGVS is in listControl. If yes, set candidate to False
        for c in listControl:
            if listAffected[0][i] in c:
                candidate = False

        #Write to file
        if candidate == True:
            listCandidate.append(listAffected[0][i].strip('\n'))

        if i%1000 == 0:
            print(str(i) + '/' + str(len(listAffected[0])) + ' completed...')

        i = i + 1

    print('Filtering completed.')
    return listCandidate

def importHGVS(listAffected, listControl):
    """
    importHGVS - imports HGVS IDs into the program for filtering
    Paramters:
    listAffected - list - list that will contain the HGVS IDs of affected individuals
    listControl - list - list that will contain the HGVS IDs of unaffected inidividuals
    Returns: none; lists passed to listAffected and listControl are populated
    """
    inputfileopen = False

    while inputfileopen == False:
        try:
            fileName = input('Please enter the name of the file listing HGVS IDs to be imported: ')
            inputFile = open(fileName, 'r')
            inputfileopen = True
        except IOError:
            print ('File not found.\n')

    cat = input('Enter "a" for HGVS IDs from affected, "c" from controls: ')

    while cat != 'a' and cat != 'c':
        cat = input('Please enter "a" for affected, and "c" for control: ')

    if cat == 'a':
        listAffected.append([])
        li = len(listAffected) - 1
        for line in inputFile:
            listAffected[li].append(line.strip('\n'))

    if cat == 'c':
        listControl.append([])
        li = len(listControl)- 1
        for line in inputFile:
            listControl[li].append(line.strip('\n'))

    print(fileName + ' successfully imported.')

def importfHGVS():
    """
    importfHGVS - imports filtered HGVS into the program
    Parameters: none
    Return: list containing the imported HGVS ids
    """

    inputfileopen = False

    while inputfileopen == False:
        try:
            fileName = input('Please enter the filename of the file listing candidate HGVS IDs: ')
            inputfile = open(fileName, 'r')
            inputfileopen = True
        except IOError:
            print('File not found.\n')

    listCandidate = []

    for line in inputfile:
        listCandidate.append(line.strip('\n'))

    return listCandidate

def mutdata(mv, HGVSID):
    """
    mutData - accepts HGVS, returns data on mutation
    Parameters: HGVSID - string - the HGVS ID of the variant to be annotated
    Return: dictionary/list containing the json data from myvariant.info
    """
    return mv.getvariant(HGVSID, fields = ['dbsnp.rsid', 'dbsnp.alleles', 'dbsnp.vartype', 'dbsnp.gene', 'clinvar.rcv.clinical_significance', 'gnomad_genome.af', 'gnomad_exome.af'])

def statCV(data):
    """
    statCV - pull ClinVar data from the JSON data structure
    Parameters: data - dictionary/list of the JSON data on the variant
    Return: the ClinVar status
    """
    pFlag = 'Benign'

    #Check if 'clinvar' exists
    if 'clinvar' in data:
        #Check if 'rcv' exists
        if 'rcv' in data['clinvar']:
            #Number of RCV entries affects parsing of the json data
            if len(data['clinvar']['rcv']) > 1:
                for y in data['clinvar']['rcv']:
                    if y['clinical_significance'] == 'Pathogenic' or y['clinical_significance'] == 'Likely pathogenic' or y['clinical_significance'] == 'Pathogenic/Likely pathogenic':
                        pFlag = 'Pathogenic'
                    elif y['clinical_significance'] == 'Uncertain significance':
                        pFlag = 'Uncertain significance'
            if len(data['clinvar']['rcv']) <= 1:
                if data['clinvar']['rcv']['clinical_significance'] == 'Pathogenic' or data['clinvar']['rcv']['clinical_significance'] == 'Likely pathogenic' or data['clinvar']['rcv']['clinical_significance'] =='Pathogenic/Likely pathogenic':
                    pFlag = 'Pathogenic'
                elif data['clinvar']['rcv']['clinical_significance'] == 'Uncertain significance':
                    pFlag = 'Uncertain significance'
            return pFlag
    else:
        return 'N/A'

def dumpCV(data):
    """
    dump_CV - pull ClinVar data from the JSON data structure
    Parameters: data - dictionary/list of the JSON data on the variant
    Return: list containing number of pathogenic, uncertain significance, and benign reports
    """
    #0 - pathogenic, 1 - uncertain significance, 2 - benign
    numReport = [0, 0, 0]
    chkFlag = False

    #Check if 'clinvar' exists
    if 'clinvar' in data:
        #Check if 'rcv' exists
        if 'rcv' in data['clinvar']:
            l = len(data['clinvar']['rcv'])
            #Number of RCV entries affects parsing of the json data
            if l > 1:
                for y in data['clinvar']['rcv']:
                    if (y['clinical_significance'] == 'Pathogenic'
                        or y['clinical_significance'] == 'Likely pathogenic'
                        or y['clinical_significance'] == 'Pathogenic/Likely pathogenic'):
                        numReport[0] = numReport[0] + 1
                    elif y['clinical_significance'] == 'Uncertain significance':
                        numReport[1] = numReport[1] + 1
                    elif (y['clinical_significance'] == 'Benign'
                          or y['clinical_significance'] == 'Likely benign'
                          or y['clinical_significance'] == 'Benign/likely benign'):
                        numReport[2] = numReport[2] + 1
                    else:
                        chkFlag = True
            if l <= 1:
                if (data['clinvar']['rcv']['clinical_significance'] == 'Pathogenic'
                    or data['clinvar']['rcv']['clinical_significance'] == 'Likely pathogenic'
                    or data['clinvar']['rcv']['clinical_significance'] =='Pathogenic/Likely pathogenic'):
                    numReport[0] = 1
                elif data['clinvar']['rcv']['clinical_significance'] == 'Uncertain significance':
                    numReport[1] = 1
                elif (data['clinvar']['rcv']['clinical_significance'] == 'Benign'
                      or data['clinvar']['rcv']['clinical_significance'] == 'Likely benign'
                      or data['clinvar']['rcv']['clinical_significance'] =='Benign/Likely benign'):
                    numReport[2] = 1
                else:
                    chkFlag = True
    #Output
    #print ('Pathogenic: ' + str(numReport[0]))
    #print ('Uncertain Significance: ' + str(numReport[1]))
    #print ('Benign: ' + str(numReport[2]))

            #if chkFlag == True:
            #   print ('Additional entry also denoted')

            return numReport

    else:
        return 'N/A'

def dumpRSID(data):
    """
    dumpRSID - pull the dbSNP RSID from the JSON data structure
    Parameters: data - dictionary/list of the JSON data on the variant
    Return: the RSID
    """
    if 'dbsnp' in data:
        return data['dbsnp']['rsid']
    else:
        return 'N/A'

def dumpvartype(data):
    """
    dumpvartype - pull the vartype from the JSON data structure
    Parameters: data - dictionary/list of the JSON data on the variant
    Return: vartype
    """
    if 'dbsnp' in data:
        return data['dbsnp']['vartype']
    else:
        return 'N/A'

def dumpgnomADG(data):
    """
    dumpgnomADG - pull the gnomAD allele frequency (inc. genome sequences) from the JSON data structure
    Parameters: data - dictionary/list of the JSON data on the variant
    Return: genome allele frequency
    """
    if 'gnomad_genome' in data:
        if 'af' in data['gnomad_genome']:
            if 'af' in data['gnomad_genome']['af']:
                return data['gnomad_genome']['af']['af']
    else:
        return 'N/A'

def dumpgnomADE(data):
    """
    dumpgnomADE - pull the gnomAD allele frequency (excl. genome sequences) from the JSON data structure
    Parameters: data - dictionary/list of the JSON data on the variant
    Return: exome allele frequency
    """
    if 'gnomad_exome' in data:
        if 'af' in data['gnomad_exome']:
            if 'af' in data['gnomad_exome']['af']:
                return data['gnomad_exome']['af']['af']
    else:
        return 'N/A'

def importanno():
    """
    importanno - from a file containing json.dumps annotation data, import it as a list of dictionaries
    Parameters: none
    Return: a list of dictionaries of annotation data
    """
    inputfileopen = False

    while inputfileopen == False:
        try:
            filename = input('Enter the name of the data file to be imported: ')
            inputfile = open(filename, 'r')
            inputfileopen = True
        except IOError:
            print('File not found.\n')

    listdata = []

    for line in inputfile:
        listdata.append(json.loads(line))

    inputfile.close()

    return listdata

def annotvep(lc):
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

    #annotlist - will store all of the data
    annotlist = []

    #Annotation process
    print('Initiating annotation...')
    while i < len(lc):

        #Retrieve 200 HGVS IDs - maximum for the batch query to the VEP REST API
        if u != (len(lc) - 1):
            data = lc[i:u]
        if u == (len(lc) - 1):
            data = lc[i:]
        pdata = str(data).replace("'", '"')

        #Retrieval of data from VEP REST API
        r = requests.post(server+ext, headers=headers, data=('{ "hgvs_notations" : ' + pdata + ' }'))

        #Error handling
        if not r.ok:
            r.raise_for_status()
            sys.exit()

        #Store as list of dictionaries
        decoded = r.json()

        #Check for missing annotations
        if len(data) != len(decoded):
            print(str(len(data) - len(decoded)) + ' annotations are missing')
            try:
                q = 0
                for hgvs in data:
                    if data.index(hgvs) < len(decoded):
                        if hgvs != decoded[q]['id']:
                            decoded.insert(data.index(hgvs), None)
                        q = q + 1

                    if data.index(hgvs) >= len(decoded):
                        decoded.append(None)
            except:
                print('Issue between ' + str(data[0]) + 'and ' + str(data[len(data)-1]))
                print(len(decoded))

        #Add these results to annot list
        for y in decoded:
            annotlist.append(y)

        print (str(u) + ' out of ' + str(len(lc)) + ' completed...')

        i = i + 200
        u = u + 200

        if u > len(lc):
            u = len(lc) - 1

    #Processing time tracker
    end = time.time()
    print('Processing time: ' + str(end - start))

    return annotlist

def annotmvi(lc):
    """
    annotmvi - from a list of HGVS IDs, retrieve all myvariant.info annotations
    Parameters:
    lc - list - list of HGVS ids
    Return: a list of the annotation data
    """
    #Import candidate mutation list
    al = []
    mv = myvariant.MyVariantInfo()

    for HGVS in lc:

        data = mv.getvariant(HGVS, fields = ['dbsnp.rsid',
                                               'dbsnp.alleles',
                                               'dbsnp.vartype',
                                               'dbsnp.gene',
                                               'clinvar.rcv.clinical_significance',
                                               'gnomad_genome.af', 'gnomad_exome.af'])
        al.append(data)

        if lc.index(HGVS)%1000 == 0:
            print (str(lc.index(HGVS)) + ' out of ' + str(len(lc)) + ' variants annotated...')

    print('All samples annotated')
    return al

def filtervariant(listanno, name = ""):
    """
    filtervariant - filters variants based on selected criteria
    Parameters: listanno - list of the annotated variants
    Returns: returns the filtered list of annotated variants
    """
    start = time.time()
    nodata = []
    nonrelevant = []
    overfivepc = []
    candidate = []
    now = datetime.datetime.now()

    if name == "":
        name = now.strftime('%y%m%d-%H%M%S')

    for data in listanno:
        nodataflag = False
        nonrelevantflag = False
        overfivepcflag = False

        if data['MScon'] == 'N/A':
            nodataflag = True
        elif (data['MScon'] == 'intron_variant'
              or data['MScon'] == 'non_coding_transcript_exon_variant'
              or data['MScon'] == 'upstream_gene_variant'
              or data['MScon'] == 'downstream_gene_variant'
              or data['MScon'] == 'synonymous_variant'
              or data['MScon'] == '5_prime_UTR_variant'
              or data['MScon'] == '3_prime_UTR_variant'
              or data['MScon'] == 'intergenic_variant'):
            nonrelevantflag = True
        elif data['gnomADG'] != 'N/A':
            if float(data['gnomADG']) >= 0.05:
                overfivepcflag = True
        elif data['gnomADE'] != 'N/A':
            if float(data['gnomADE']) >= 0.05:
                overfivepcflag = True

        if nodataflag == True:
            nodata.append(data)
        elif nonrelevantflag == True:
            nonrelevant.append(data)
        elif overfivepcflag == True:
            overfivepc.append(data)
        else:
            candidate.append(data)
        if listanno.index(data)%1000 == 0:
            print (str(listanno.index(data)) + ' out of ' + str(len(listanno)) + ' filtered...')

    #Write to nodatavariant.txt
    outputfilendv = open('nodatavariant_' + name + '.txt', 'w')
    outputfilendv.write('HGVS' + '\t'
                       + 'RSID' + '\t'
                       + 'vartype' + '\t'
                       + 'gnomADG' + '\t'
                       + 'gnomADE' + '\t'
                       + 'ClinVar' + '\t'
                       + 'MScon' + '\n')
    for data in nodata:
        outputfilendv.write(data['_id'] + '\t'
                         + data['rsid'] + '\t'
                         + data['vartype'] + '\t'
                         + str(data['gnomADG']) + '\t'
                         + str(data['gnomADE']) + '\t'
                         + str(data['ClinVar']) + '\t'
                         + data['MScon'] + '\n')
    outputfilendv.close()

    #Write to nonrelevant.txt
    outputfilenr = open('nonrelevantvariant_' + name + '.txt', 'w')
    outputfilenr.write('HGVS' + '\t'
                       + 'RSID' + '\t'
                       + 'vartype' + '\t'
                       + 'gnomADG' + '\t'
                       + 'gnomADE' + '\t'
                       + 'ClinVar' + '\t'
                       + 'MScon' + '\n')
    for data in nonrelevant:
        outputfilenr.write(data['_id'] + '\t'
                         + data['rsid'] + '\t'
                         + data['vartype'] + '\t'
                         + str(data['gnomADG']) + '\t'
                         + str(data['gnomADE']) + '\t'
                         + str(data['ClinVar']) + '\t'
                         + data['MScon'] + '\n')
    outputfilenr.close()

    #Write to overfivepc.txt
    outputfileofp = open('overfivepcvariant_' + name + '.txt', 'w')
    outputfileofp.write('HGVS' + '\t'
                       + 'RSID' + '\t'
                       + 'vartype' + '\t'
                       + 'gnomADG' + '\t'
                       + 'gnomADE' + '\t'
                       + 'ClinVar' + '\t'
                       + 'MScon' + '\n')
    for data in overfivepc:
        outputfileofp.write(data['_id'] + '\t'
                         + data['rsid'] + '\t'
                         + data['vartype'] + '\t'
                         + str(data['gnomADG']) + '\t'
                         + str(data['gnomADE']) + '\t'
                         + str(data['ClinVar']) + '\t'
                         + data['MScon'] + '\n')
    outputfileofp.close()

    #Write to candidatevariants.txt
    outputfile = open('candidatevariants_' + name + '.txt', 'w')
    outputfile.write('HGVS' + '\t'
                       + 'RSID' + '\t'
                       + 'vartype' + '\t'
                       + 'gnomADG' + '\t'
                       + 'gnomADE' + '\t'
                       + 'ClinVar' + '\t'
                       + 'MScon' + '\n')
    for data in candidate:
        outputfile.write(data['_id'] + '\t'
                         + data['rsid'] + '\t'
                         + data['vartype'] + '\t'
                         + str(data['gnomADG']) + '\t'
                         + str(data['gnomADE']) + '\t'
                         + str(data['ClinVar']) + '\t'
                         + data['MScon'] + '\n')
    outputfile.close()

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
                        + '# Greater than 5% AF: '
                        + str(len(overfivepc))
                        + '\n'
                        + '# Candidates: '
                        + str(len(candidate))
                        + '\n'
                        + 'Processing time: '
                        + str (end - start)
                        + ' seconds'
                        )

    return candidate

def importmut():
    """
    importmut - reads file containing annotations generated by this program to repopulate the data structure
    Parameters: none
    Returns: a list of dictionaries containing all of the annotation data
    """
    listanno = []
    inputfileopen = False

    while inputfileopen == False:
        try:
            filename = input ('Please enter the name of the file containing the annotated mutations: ')
            inputfile = open(filename, 'r')
            inputfileopen = True
        except IOError:
            print ('File not found.')

    #Skip header
    next(inputfile)

    for line in inputfile:
        data = line.split('\t')
        listanno.append(dict({'_id':data[0],
                           'rsid':data[1],
                           'vartype':data[2],
                           'gnomADG':data[3],
                           'gnomADE':data[4],
                           'ClinVar':data[5],
                           'MScon':data[6].strip('\n')
                           }))

    inputfile.close()

    print('Import completed.')

    return listanno

def combineanno(listmvi, listvep):
    """

    combineanno - combines annotations from MVI and VEP together into one list of dictionaries
    Parameters:
    listmvi - list - contains a list of dictionaries of MVI information
    listvep - list - contains a list of dictionaries of VEP information
    Returns: a list containing both sets of information

    """
    al = []
    lc = importfHGVS()
    i = 0

    while i < len(lc):

        if str(type(listvep[i])) != "<class 'NoneType'>":
            consequence = listvep[i]['most_severe_consequence']
        else:
            consequence = 'N/A'

        if str(type(listmvi[i])) != "<class 'NoneType'>":
            al.append(dict({'_id':listmvi[i]['_id'],
                            'rsid':dumpRSID(listmvi[i]),
                            'vartype':dumpvartype(listmvi[i]),
                            'gnomADG':dumpgnomADG(listmvi[i]),
                            'gnomADE':dumpgnomADE(listmvi[i]),
                            'ClinVar':dumpCV(listmvi[i]),
                            'MScon':consequence
                            }))

        if str(type(listmvi[i])) == "<class 'NoneType'>":
            al.append(dict({'_id':lc[i],
                            'rsid':'N/A',
                            'vartype':'N/A',
                            'gnomADG':'N/A',
                            'gnomADE':'N/A',
                            'ClinVar':'N/A',
                            'MScon':consequence
                            }))

        i = i + 1

    return al

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
                     + 'MSCon' + '\n')

    #Write data
    for i in listanno:
        outputfile.write(i['_id'] + '\t'
                         + i['rsid'] + '\t'
                         + i['vartype'] + '\t'
                         + str(i['gnomADG']) + '\t'
                         + str(i['gnomADE']) + '\t'
                         + str(i['ClinVar']) + '\t'
                         + i['MScon'])
        if listanno.index(i) != (len(listanno)-1):
            outputfile.write('\n')

        #Progress Indicator
        if listanno.index(i)%1000 == 0:
            print (str(listanno.index(i)) + ' out of ' + str(len(listanno)) + ' written...')

    print('All annotations written to file')

    outputfile.close()

def exportanno(listanno, filename):
    """
    exportanno - exports raw annotation data
    Parameters:
    listanno - list - list of dictionaries containing raw JSON data
    filename - string - name of the file to output to
    Returns: none; outputs to file with filename
    """
    outputfile = open(filename, 'w')
    for y in listanno:
        outputfile.write(json.dumps(y))
        if listanno.index(y) != (len(listanno)-1):
            outputfile.write('\n')
    outputfile.close()

def retrievevep(lc, listvep):
    al = []
    i = 0
    vi = 0

    while i < len(lc):
        found = False

        if i%10 == 0:
            print(str(i) + ' out of ' + str(len(lc)) + ' completed....')

        while (found == False) and vi < len(listvep):
            if lc[i]['_id'] == listvep[vi]['id']:
                found = True
                al.append(listvep[vi])
            vi = vi + 1

        if vi > len(listvep):
            if i < len(listvep):
               vi = i
            else:
               vi = 19000

        i = i + 1

    print('Samples completed.')
    return al

def transcriptids(listanno):
    listall = []
    listrelevant = []

    for variant in listanno:
        #Lists of transcripts for the variants
        list_relevant_transcript = []
        list_transcript = []

        print ('Sample ' + str(listanno.index(variant)))
        #Check if data has been pulled
        if 'transcript_consequences' in variant:
            #Go through each transcript
            for transcript in variant['transcript_consequences']:
                #Flag for relevance
                append_relevant = False
                #Go through each consequence term
                for terms in transcript['consequence_terms']:
                    #Check if term is relevant
                    if (terms != 'intron_variant'
                        and terms != 'non_coding_transcript_exon_variant'
                        and terms != 'upstream_gene_variant'
                        and terms != 'downstream_gene_variant'
                        and terms != 'synonymous_variant'
                        and terms != '5_prime_UTR_variant'
                        and terms != '3_prime_UTR_variant'
                        and terms != 'intergenic_variant'
                        and terms != 'non_coding_transcript_variant'):

                        append_relevant = True

                #Final check if the protein is coding
                if (transcript['biotype']) != 'protein_coding':
                    append_relevant = False

                #Verify information is present
                protein_start = 'None'
                protein_end = 'None'
                amino_acids = 'None'

                if 'protein_start' in transcript:
                    protein_start = transcript['protein_start']
                if 'protein_end' in transcript:
                    protein_end = transcript['protein_end']
                if 'amino_acids' in transcript:
                    amino_acids = transcript['amino_acids']

                #Relevant
                if append_relevant == True:
                    list_relevant_transcript.append(dict({'transcript_id':transcript['transcript_id'],
                                                         'consequence_terms':transcript['consequence_terms'],
                                                         'biotype':transcript['biotype'],
                                                         'protein_start':protein_start,
                                                         'protein_end':protein_end,
                                                         'amino_acids':amino_acids
                                                         }))
                list_transcript.append(dict({'transcript_id':transcript['transcript_id'],
                                            'consequence_terms':transcript['consequence_terms'],
                                            'biotype':transcript['biotype'],
                                            'protein_start':protein_start,
                                            'protein_end':protein_end,
                                            'amino_acids':amino_acids
                                            }))

        else:
            list_relevant_transcript.append('None')
            list_transcript.append('None')

        listall.append(dict({'id':variant['id'],
                             'transcripts':list_transcript,
                             }))
        listrelevant.append(dict({'id':variant['id'],
                                  'transcripts':list_relevant_transcript,
                                  }))

    exportanno(listall, 'variants_all_transcripts.txt')
    exportanno(listrelevant, 'variants_relevant_transcript.txt')

def vcftoHGVS():
    HGVS = []
    filename = input('Please enter filename : ')
    HGVS = list(myvariant.get_hgvs_from_vcf(filename))

    index = 0

    #Accept name input from user
    name = input('Please enter an identifier for your file: ')
    outputfile = open('HGVS_' + name + '.txt', 'w')

    print('Writing to file...')
    #Write to file
    for ids in HGVS:
        outputfile.write(ids + '\n')
    print('Done.')

    outputfile.close()
