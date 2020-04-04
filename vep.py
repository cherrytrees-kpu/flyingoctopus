"""

vep.py - code for processing data from Ensembl's Variant Effect Predictor

"""
#Imported modules
import time
import datetime
import requests

#Function definitions
def annotate(listHGVS):
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
    parameters = {"canonical":"1", "uniprot":"1"}

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
        r = requests.post(server+ext, headers=headers, data=('{ "hgvs_notations" : ' + p_rangeHGVS + ' }'), params = parameters)

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

def dumpensembltranscripts(anno_vep):
    list_return = []

    #Check if data is actually present
    if anno_vep is not None and 'transcript_consequences' in anno_vep:
        #Identify number of genes that the variant impacts
        #Create list of genes that variant impacts
        for transcript in anno_vep['transcript_consequences']:
            gene = {'gene_id':transcript['gene_id'],
                    'gene_symbol':transcript['gene_symbol'],
                    'transcripts':[]
                    }
            if gene not in list_return:
                list_return.append(gene)

        #Go through all of the transcript the variant impacts
        for transcript in anno_vep['transcript_consequences']:
            #If the transcript is relevant, add it to appropriate gene
            relevantflag = False
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
            if relevantflag is True:
                for gene in list_return:
                    if transcript['gene_id'] == gene['gene_id']:
                        gene['transcripts'].append(transcript)
        #Remove genes with no relevant transcripts
        i = 0
        while i < len(list_return):
            if len(list_return[i]['transcripts']) == 0:
                del list_return[i]
                i = i - 1
            i = i + 1

        #Return the list
        return list_return
    else:
        return None

def ensemblsequence(transcriptid, seq):
    server = 'https://grch37.rest.ensembl.org'
    ext = '/sequence/id/' + transcriptid + '?'

    xheaders={'Content-Type':'application/json'}
    xparams = {'type':seq,'object_type':'transcript'}

    r = requests.get(server+ext, headers=xheaders, params=xparams)

    if r.status_code == 200:
        return r.json()['seq']
    else:
        print(seq + ' for ' + transcriptid + ' could not be found.')
        return None
