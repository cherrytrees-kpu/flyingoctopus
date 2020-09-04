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
