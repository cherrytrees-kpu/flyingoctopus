"""

mvi.py - code for processing data from myvariant.info

"""
#Import modules
import myvariant

#Function definitions
def annotate(listHGVS):
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

def getmvi(idHGVS, mv = myvariant.MyVariantInfo()):
    anno = mv.getvariant(idHGVS, fields = ['dbsnp.rsid',
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
                                        ])
    return anno

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
