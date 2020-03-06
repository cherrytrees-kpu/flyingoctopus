"""

hpa.py - code for processing data from the Human Protein Atlas

"""
#Import modules
import requests

#Function definitions
def annotate(geneid):
    """
    annothpa - pulls expression data from the Human Protein Atlas
    Parameters: geneid - Ensembl gene id
    Returns: returns expression data
    """
    #expression - list that holds the return data from

    server = "http://www.proteinatlas.org/"

    r = requests.get(server + geneid + '.json')

    #Check if there was data returned
    if r.status_code == requests.codes.ok:
        return r.json()
    else:
        return None
