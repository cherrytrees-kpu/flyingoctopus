"""

vep.py - code for processing data from Ensembl's Variant Effect Predictor

"""
#Imported modules
import time
from ratelimit import limits, sleep_and_retry
import datetime
import requests

#Function definitions
@sleep_and_retry
@limits(calls=55000, period=3600)
def getvep(idHGVS):
    #VEP REST API
    server = "http://grch37.rest.ensembl.org"
    ext = "/vep/human/hgvs/"
    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
    parameters = {"canonical":"1", "uniprot":"1"}

    anno = None

    #Get the annotation
    try:
        r = requests.get(server+ext+idHGVS, headers=headers, params = parameters)
    except ConnectionError:
        print('Connection failed! Data for variant not retrieved.')
        print('Trying again after 10 seconds.....')
        time.sleep(10)
        try:
            r = requests.get(server+ext+idHGVS, headers=headers, params = parameters)
        except ConnectionError:
            print('Connection failed again! Data for variant not retrieved, try again later.')
    except requests.exceptions.HTTPError:
        print('No data available for '+ idHGVS +'.')
    else:
        anno = r.json()
    finally:
        return anno

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

def parsevepanno(data_variant, vepanno):
    def getconsequence(vepanno):
        """
        getconsequence - pull the most severe consequence from the JSON data structure
        Parameters: anno_vep - dictionary/list of the JSON data on the variant
        Return: most severe consequence
        """
        return vepanno['most_severe_consequence']
    def getgenes(vepanno):
        #Determine if the gene does affect transcripts
        if 'transcript_consequences' in vepanno:
            set_genes = set()
            for transcript in vepanno['transcript_consequences']:
                set_genes.add(transcript['gene_id']+'/'+transcript['gene_symbol'])
            return list(set_genes)
        else:
            print('Does not affect gene sequence.')

    data_variant['MScon'] = dumpconsequence(vepanno)
    listgenes = getgenes(vepanno)
    data_variant['genes'] = []
    for gene in listgenes:
        data_variant['genes'].append({'gene':gene.split('/')[1], 'gene_id':gene.split('/')[0]})
