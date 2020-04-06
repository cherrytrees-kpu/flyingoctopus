"""

hpa.py - code for processing data from the Human Protein Atlas

"""
#Import modules
import requests

#Function definitions
def gethpa(geneid):
    """
    annothpa - pulls expression data from the Human Protein Atlas
    Parameters: geneid - Ensembl gene id
    Returns: returns expression data
    """

    server = "http://www.proteinatlas.org/"

    anno = None

    #Get the annotation
    try:
        r = requests.get(server + geneid + '.json')
    except ConnectionError:
        print('Connection failed! Data for gene not retrieved.')
        print('Trying again after 10 seconds.....')
        time.sleep(10)
        try:
            r = requests.get(server + geneid + '.json')
        except ConnectionError:
            print('Connection failed again! Data for gene not retrieved, try again later.')
    except requests.exceptions.HTTPError:
        print('No HPA data available for '+ geneid +'.')
    else:
        anno = r.json()
    finally:
        return anno

def annotate(annotated_variants, getfromfile = False):
    """
    annotate - pulls expression data from the Human Protein Atlas
    Parameters: listvariants - list of variants
    Returns: returns with data
    """
    if getfromfile == False:
        print('yes')
