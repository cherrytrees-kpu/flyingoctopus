"""

uniprot.py - code for Uniprot data processing and handling

"""
import requests
from ratelimit import limits, sleep_and_retry
import xml.etree.ElementTree as ET

@sleep_and_retry
@limits(calls=1, period=10)
def retrievexml(id_uniprot):
    #Server
    server = 'http://uniprot.org/uniprot/'
    #Pulling data
    r = requests.get(server + id_uniprot + '.xml')
    #Checking if all is good
    if r.status_code == requests.codes.ok:
        print('Data successfully accessed.')
    else:
        print('Data not retrieved.')
    return(r.text)

def uniprotparse(xml):
    #Generating XML tree
    tree = ET.fromstring(xml)
    features = []
    #Go through all of the features
    for element in tree[0].iter('{http://uniprot.org/uniprot}feature'):
        feature = element.attrib
        x = element.findall('.*')
        for e in x:
            name = (e.tag).split('}')[1]
            feature[name] = e.text
            if name == 'location':
                y = e.findall('.*')
                for pos in y:
                    name = (pos.tag).split('}')[1]
                    feature[name] = pos.get('position')
        #Delete the empty location tag
        del feature['location']
        features.append(feature)

    #Pull the function of the protein
    function = ''
    for element in tree[0].iter('{http://uniprot.org/uniprot}comment'):
        comment = element.attrib
        if comment['type'] == 'function':
            function = element.find('{http://uniprot.org/uniprot}text').text

    #Pull protein name
    try:
        protein_name = tree[0].find('{http://uniprot.org/uniprot}protein').find('{http://uniprot.org/uniprot}recommendedName').find('{http://uniprot.org/uniprot}fullName').text
    except AttributeError:
        print ('No protein name')
        protein_name = None
    #Pull dataset
    dataset = tree[0].attrib.get('dataset')
    data = {'protein_name':protein_name, 'function':function,'features':features, 'dataset':dataset}
    return data

def uniprotacc(geneid):
    #server
    server = 'https://www.uniprot.org/uploadlists'
    map_parameters = {
    'from':'ENSEMBL_ID',
    'to':'ACC',
    'format':'tab',
    'query':geneid,
    }
    #Get the Uniprot IDs
    r = requests.get(server, params = map_parameters)
    #Parse data
    if r.status_code == requests.codes.ok:
        data = ((r.text).strip('\n')).split('\n')
        del data[0]
        list_UniP = []
        for entry in data:
            list_UniP.append(entry.split('\t')[1])
        print ('IDs are mapped.')
        return list_UniP
    else:
        print('ID did not map.')

def annotate(gene_id):
    #Variables
    list_data = []
    list_unip = uniprotacc(gene_id)
    print(gene_id)
    #Get Uniprot data
    for unip in list_unip:
        xml = retrievexml(unip)
        data = uniprotparse(xml)
        data['acc'] = unip
        data['gene_id'] = 'gene_id'
        if data['dataset'] == 'Swiss-Prot':
            print('Added to data list')
            list_data.append(data)
        else:
            print('Not in Swiss-Prot')
    if len(list_data) == 1:
        return list_data[0]
    elif len(list_data) == 0:
        return None
    else:
        print('Multiple entries found for ' + str(gene_id))
