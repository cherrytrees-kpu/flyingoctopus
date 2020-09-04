"""

mvi.py - code for processing data from myvariant.info

"""
#Import modules
import myvariant
import requests

#Function definitions
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

#Function definitions
def annotateHPA(geneid):
    """
    annothpa - pulls expression data from the Human Protein Atlas
    Parameters: geneid - Ensembl gene id
    Returns: returns expression data
    """
    #expression - list that holds the return data from

    server = "http://www.proteinatlas.org/"

    try:
        r = requests.get(server + geneid + '.json', timeout=20)
    except requests.exceptions.ConnectionError:
        print('ConnectionError')
        try:
            r = requests.get(server + geneid + '.json', timeout=20)
        except requests.exceptions.ConnectionError:
            print('Connection Error again')
            return "connectionError"
        else: 
            return r.json()
    except requests.exceptions.TooManyRedirects:
        print('Too Many Redirects')
        return 'tooManyRedirects'
    else:
        return r.json()

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
