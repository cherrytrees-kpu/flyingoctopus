"""
VariantAnnotation.py - handles variant annotation tasks

Author: Michael Ke

"""
import json
import sys
import pathlib
import shutil
import requests
import myvariant
import time
import va
import datetime
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def fullroutine():
    #Get unique job ID name
    now = datetime.datetime.now()
    timelabel = now.strftime('%y%m%d-%H%M%S')

    #Import candidate variants
    lCandidate = va.importfHGVS()

    #MVI Annotation
    startmvi = time.time()
    lmvi = va.annotmvi(lCandidate)
    endmvi = time.time()
    print('Total MVI annotation run time: ' + str(endmvi - startmvi) + '\n')
    #Export MVI
    mvifname = 'mvi_annotation_' + timelabel + '.txt'
    va.exportanno(lmvi, mvifname)
    #VEP annotation
    startvep = time.time()
    lvep = va.annotvep(lCandidate)
    endvep = time.time()
    print ('Total VEP annotation run time: '+  str(endvep - startvep) + '\n')
    #Export VEP
    vepfname = 'vep_annotation_' + timelabel + '.txt'
    va.exportanno(lvep, vepfname)
    #Export combined annotations
    lanno = va.combineanno(lmvi, lvep)
    annofname = 'annotated_mutations_' + timelabel + '.txt'
    va.writeanno(lanno, annofname)
    #Filtervariant
    va.filtervariant(lanno, timelabel)
    print('Complete')

def filternodata(listanno):
    list_nodata = []
    list_candidate = []

    for anno in listanno:
        nodata = va.checknodata(anno)
        if nodata is True:
            list_nodata.append(anno)
        elif nodata is False:
            list_candidate.append(anno)

    return list_nodata, list_candidate

def filterfreq(listanno):
    list_overfreq = []
    list_candidate = []

    for anno in listanno:
        overfreq = va.checkfreq(anno)
        if overfreq is True:
            list_overfreq.append(anno)
        elif overfreq is False:
            list_candidate.append(anno)

    return list_overfreq, list_candidate

def filtercons(listanno):
    list_irrelevant = []
    list_candidate = []

    for anno in listanno:
        irrelevant = va.checkcons(anno)
        if irrelevant is True:
            list_irrelevant.append(anno)
        elif irrelevant is False:
            list_candidate.append(anno)

    return list_irrelevant, list_candidate

def filterexpression(listanno):
    list_hpa = []
    list_notbrainexpressed = []
    list_candidate = []

    #Get annotation from Human Protein Atlas
    for anno in listanno:
        if anno['genes'] is not None:
            for gene in anno['genes']:
                anno_hpa = va.annothpa(gene['gene_id'])
                list_hpa.append(anno_hpa)
                if anno_hpa is not None:
                    gene['RNAbrd'] = anno_hpa['RNA brain regional distribution']
                else:
                    gene['RNAbrd'] = None

    #Filtering
    for anno in listanno:
        notbrainexpressflag = True
        if anno['genes'] is not None:
            for gene in anno['genes']:
                if gene['RNAbrd'] != 'Not detected':
                    notbrainexpressflag = False

        if notbrainexpressflag is True:
            list_notbrainexpressed.append(anno)
        else:
            list_candidate.append(anno)

    #Export and return
    va.exportanno(list_hpa, 'hpa_annotations.txt')
    return list_notbrainexpressed, list_candidate

def getsequences(listanno, basepath):
    def substitution(transcript, pos, allele, strand):
        if strand == -1:
            allele = str(Seq(allele, IUPAC.unambiguous_dna).complement())
        return transcript[:pos-1] + allele + transcript[pos:]

    def insertion(transcript, start, end, insert, strand):
        if strand == -1:
            insert = str(Seq(insert, IUPAC.unambiguous_dna).complement())
        return transcript[:start] + insert + transcript[end:]

    def deletion(transcript, start, end):
        return transcript[:start-1] + transcript[end:]

    def translation (sequence):
        if 'N' not in sequence:
            dna = Seq(sequence, IUPAC.unambiguous_dna)
            mrna = dna.transcribe()
            protein = mrna.translate()
            return str(protein)
        else:
            print('Ambiguous nucleotide present, no protein')

    #Do for every variant
    for anno in listanno:
        #Holds stats
        list_stats = []
        #Make a folder for variant
        print ('Working on ' + anno['_id'] + '...')
        variantpath = basepath.joinpath(va.converthgvs(anno['_id']))
        print(variantpath)
        variantpath.mkdir(exist_ok=True)
        #Do for every gene
        for gene in anno['genes']:
            #Stats to track
            list_transcripts = []
            list_notranscripts = []
            #Make folder for gene
            genepath = variantpath.joinpath(gene['gene_id'])
            genepath.mkdir(exist_ok=True)
            #Do for every transcript
            for transcript in gene['transcripts']:
                list_transcripts.append(transcript['transcript_id'])
                print('Pulling transcripts for ' + transcript['transcript_id'])
                cdna = va.ensemblsequence(transcript['transcript_id'], 'cdna')
                cds = va.ensemblsequence(transcript['transcript_id'], 'cds')
                #Write cdna and cds
                va.outputsequence(cdna, transcript['transcript_id'], 'cdna', genepath)
                va.outputsequence(cds, transcript['transcript_id'], 'cds', genepath)
                #Get rid of splice_region_variant, TF_binding_site, regulatory_region_variant
                if cds is not None:
                    get_protein = True
                else:
                    get_protein = False
                    list_notranscripts.append(transcript['transcript_id'])
                for term in transcript['consequence_terms']:
                    if (term == 'splice_region_variant'
                        or term == 'TF_binding_site'
                        or term == 'regulatory_region_variant'):
                        get_protein = False
                #Get protein sequences
                if (get_protein is True) and (anno['vartype'] == 'snv'):
                    vcds = substitution(cds,
                                        transcript['cds_start'],
                                        transcript['variant_allele'],
                                        transcript['strand'],
                                        )
                    vprotein = translation(vcds)
                elif (get_protein is True) and (anno['vartype'] == 'del'):
                    vcds = deletion(cds,
                                    transcript['cds_start'],
                                    transcript['cds_end'],
                                    )
                    vprotein = translation(vcds)
                elif (get_protein is True) and (anno['vartype'] == 'ins' or anno['vartype' == 'delins']):
                    vcds = insertion(cds,
                                    transcript['cds_start'],
                                    transcript['cds_end'],
                                    transcript['variant_allele'],
                                    transcript['strand']
                                    )
                    vprotein = translation(vcds)
                #Write proteins
                if get_protein is True:
                    protein = translation(cds)
                    va.outputsequence(protein, transcript['transcript_id'], 'protein', genepath)
                    va.utputsequence(vprotein, transcript['transcript_id'], 'vprotein', genepath)
            #Add to the stats document
            list_stats.append(dict({'gene':gene['gene_id'],
                                    'transcripts':list_transcripts,
                                    'no_transcripts':list_notranscripts,
                                    }))
        #Write the summary file for the variant
        summary_file = open(variantpath.as_posix()
                            + '/'
                            + va.converthgvs(anno['_id'])
                            + '_summary.txt'
                            ,
                            'w',
                            )
        for item in list_stats:
            summary_file.write(item['gene']
                                + '\n'
                                + 'All transcripts:'
                                + str(item['transcripts'])
                                + '\n'
                                + 'No sequence:'
                                + str(item['no_transcripts'])
                                + '\n'
                                )

        summary_file.close()

##### Program Start ##################################################################
#GLOBAL VARIABLES
EXIT_PROGRAM = False
#LIST_ANNO - holds all annotations
LIST_ANNO = []
#LIST_VEP - holds VEP annotations
LIST_VEP = []
#LIST_MVI - holds myvariant.info annotations
LIST_MVI = []
#LIST_CANDIDATE = holds the candidate annotations
LIST_CANDIDATE = []

BASEPATH = pathlib.Path.cwd()

while EXIT_PROGRAM == False:
    #Display menu
    print('1) Import HGVS from .vcf files')
    print('2) Generate filtered HGVS file')
    print('3) Annotation')
    print('4) Import data')
    print('5) Export data')
    print('6) Process data')
    print('7) Get Ensembl CDS sequences')
    print('8) Exit (modules still loaded)')

    #Accept input from user
    option = int(input('Select option: '))
    #Input check
    while (option != 1
           and option != 2
           and option != 3
           and option != 4
           and option != 5
           and option != 6
           and option != 7
           and option != 8
           ):
        option = int(input('Invalid selection; please select one of the options: '))
    print('')

    #1) Import HGVS from .vcf files
    if option == 1:

        filename = input('Please enter the .vcf file to extract HGVS from: ')
        va.vcftoHGVS(filename)

    #2) Generate filtered HGVS file
    if option == 2:
        #Variables
        start = time.time()
        listaffected = []
        listcontrol = []
        numfiles = int(input('How many files individuals would you like to use: '))
        i = 0

        #Import the number of files that will be analyzed
        while i < numfiles:
            #Open the HGVS ID file
            try:
                filename = input('Please enter the name of the file listing HGVS IDs to be imported: ')
                inputfile = open(filename, 'r')
            except IOError:
                print('File not found.')
                filename = input('Re-enter filename: ')

            #Ask user if the HGVS ID file is from an affected or control
            cat = input('Enter "a" for HGVS IDs from affected, "c" from controls: ')
            while cat != 'a' and cat != 'c':
                cat = input('Please enter "a" for affected, and "c" for control: ')

            #For affected:
            if cat == 'a':
                listaffected.append([])
                li = len(listaffected) - 1
                for idHGVS in inputfile:
                    listaffected[li].append(idHGVS.strip('\n'))
            #For control:
            if cat == 'c':
                listcontrol.append([])
                li = len(listcontrol)- 1
                for line in inputfile:
                    listcontrol[li].append(idHGVS.strip('\n'))
            print(filename + ' successfully imported.')
            i = i + 1

        #Perform filtering
        listfiltered= va.filteraffected(listaffected, listcontrol)
        list_gl = []

        #Remove any GL variants from the filtered list
        j = 0
        while j < len(listfiltered):
            if 'chrGL' in listfiltered[j]:
                list_gl.append(listfiltered[j])
                del listfiltered[j]
                j = j - 1
            j = j + 1

        va.outputHGVS(lCandidate, "candidate")
        va.outputHGVS(lglvariant, "GLvariant")
        end = time.time()
        print('Total run time: ' + str(end - start) + '\n')

    #3) Annotation
    elif option == 3:
        exitannotation = False
        while exitannotation == False:
            print ('1) myvariant.info annotation')
            print ('2) VEP annotation')
            print ('3) Return to main menu')

            optionannotation = int(input('Select option: '))

            while (optionannotation != 1
                   and optionannotation != 2
                   and optionannotation != 3):
                optionannotation = int(input('Invalid selection; please select one of the options: '))

            if optionannotation == 1:
                #Open HGVS ID file that will be annotated
                filename = input('Please enter the name of the file listing HGVS IDs to be annotated: ')
                listHGVS = va.importHGVS(filename)
                #Perform MVI annotation
                start = time.time()
                LIST_MVI = va.annotmvi(listHGVS)
                end = time.time()
                print('Total run time: ' + str(end - start) + '\n')

            elif optionannotation == 2:
                #Open HGVS ID file that will be annotated
                filename = input('Please enter the name of the file listing HGVS IDs to be annotated: ')
                listHGVS = va.importHGVS(filename)
                #Perform VEP annotation
                start = time.time()
                LIST_VEP = va.annotvep(listHGVS)
                end = time.time()
                print('Total run time: ' + str(end - start) + '\n')

            elif optionannotation == 3:
                exitannotation = True
                print('')

    elif option == 4:
        exitimport = False
        while exitimport == False:
            print ('1) Import myvariant.info annotation data')
            print ('2) Import VEP annotation data')
            print ('3) Import program-generated data')
            print ('4) Return to main menu')

            optionimport = int(input('Select option: '))

            while (optionimport != 1
                   and optionimport !=2
                   and optionimport !=3
                   and optionimport !=4):
                optionimport = int(input('Invalid selection; please select one of the options: '))

            if optionimport == 1:
                #Open the HGVS ID file
                filename = input('Please enter the name of the file being imported : ')
                LIST_MVI = va.importanno(filename)
                print('Data imported.' + '\n')
            elif optionimport == 2:
                filename = input('Please enter the name of the file being imported : ')
                LIST_VEP = va.importanno(filename)
                print('Data imported.' + '\n')
            elif optionimport == 3:
                filename = input('Please enter the name of the file being imported : ')
                LIST_ANNO = va.importanno(filename)
                print('Data imported.' + '\n')
            elif optionimport == 4:
                exitimport = True
                print('')

    elif option == 5:
        exitexport = False
        while exitexport == False:
            print('1) Export myvariant.info annotation data')
            print('2) Export VEP annotation data')
            print('3) Export combined data')
            print('4) Return to main menu')

            optionexport = int(input('Select option: '))

            while (optionexport != 1
                   and optionexport != 2
                   and optionexport != 3
                   and optionexport != 4):
                optionexport = int(input('Invalid selection; please select one of the options: '))

            if optionexport == 1:
                va.exportanno(LIST_MVI, 'myvariantannotation.txt')
            elif optionexport == 2:
                va.exportanno(LIST_VEP, 'vepannotation.txt')
            elif optionexport == 3:
                filename = input('Please enter the name of the file being imported : ')
                listHGVS = va.importHGVS(filename)
                LIST_ANNO = va.combineanno(LIST_MVI, LIST_VEP, listHGVS)
                va.exportanno(LIST_ANNO, 'annotated_mutations.txt')
            elif optionexport == 4:
                exitexport = True
                print('')

    elif option == 6:
        #Lists
        list_nodata = []
        list_irrelevant = []
        list_highfreq = []
        list_notexpressedbrain = []
        list_candidate = []
        list_intermediate = []

        #Job Identification
        now = datetime.datetime.now()
        time = now.strftime('%y%m%d-%H%M%S')

        #File-handling variables
        basepath = pathlib.Path.cwd()
        file_created = False
        while file_created is False:
            file_ID = input('Enter a job ID: ')
            newpath = basepath.joinpath(file_ID+'_'+time)
            try:
                newpath.mkdir()
            except FileExistsError as fee_error:
                print(fee_error)
                print('Please try another ID: ')
            else:
                file_created = True

        #Filenames
        filename_nodata = 'nodata_' + str(file_ID) + '_' + time + '.txt'
        filename_irrelevant = 'nonrelevant_' + str(file_ID) + '_' + time + '.txt'
        filename_highfreq = 'overfreqpc_' + str(file_ID) + '_' + time + '.txt'
        filename_notexpressedbrain = 'notbrainexpress_' + str(file_ID) + '_' + time + '.txt'
        filename_candidate = 'candidatevariants_' + str(file_ID) + '_' + time + '.txt'
        filename_summary = 'summary_' + str(file_ID) + '_' + time + '.txt'

        #Conduct filtering
        #No data
        print('Filtering out no data...')
        list_nodata, list_intermediate = filternodata(LIST_ANNO)
        print('No data variants filtered out.')
        #Consequence filter
        print('Filtering out irrelevant variants...')
        list_irrelevant, list_intermediate = filtercons(list_intermediate)
        print('Irrelevant variants filtered out.')
        #Frequency filter
        print('Filtering out high frequency variants...')
        list_highfreq, list_intermediate = filterfreq(list_intermediate)
        print('High frequency variants filtered out.')
        #Expression filter
        print('Filtering out variants not expressed in the brain...')
        list_notexpressedbrain, list_candidate = filterexpression(list_intermediate)
        print('Variants not expressed in brain filtered out.')

        LIST_CANDIDATE = list_candidate

        #Export all of the files
        va.exportanno(list_nodata, newpath.as_posix() + '/' + filename_nodata)
        va.exportanno(list_irrelevant, newpath.as_posix() + '/' + filename_irrelevant)
        va.exportanno(list_highfreq, newpath.as_posix() + '/' + filename_highfreq)
        va.exportanno(list_notexpressedbrain, newpath.as_posix() + '/' + filename_notexpressedbrain)
        va.exportanno(list_candidate, newpath.as_posix() + '/' + filename_candidate)

        #Export candidate by snp type
        list_snv = []
        list_del = []
        list_ins = []
        for candidate in list_candidate:
            if candidate['vartype'] == 'snv':
                print('snp')
                list_snv.append(candidate)
            elif candidate['vartype'] == 'del':
                print ('del')
                list_del.append(candidate)
            elif candidate['vartype'] == 'ins' or candidate['vartype'] == 'delins':
                print ('delins')
                list_ins.append(candidate)

        va.exportanno(list_snv,
                        newpath.as_posix()
                        + '/'
                        + 'candidatesnv_'
                        + str(file_ID)
                        + '_'
                        + time
                        + '.txt'
                        )
        va.exportanno(list_del,
                        newpath.as_posix()
                        + '/'
                        + 'candidatedel_'
                        + str(file_ID)
                        + '_'
                        + time
                        + '.txt'
                        )
        va.exportanno(list_ins,
                        newpath.as_posix()
                        + '/'
                        + 'candidateins_'
                        + str(file_ID)
                        + '_'
                        + time
                        + '.txt'
                        )
        #Summary
        summaryfile = open(newpath.as_posix() + '/' + filename_summary, 'w')

        summaryfile.write('Summary report of analysis'
                        + '\n'
                        + 'ID: '
                        + file_ID
                        + '\n'
                        + 'Date performed: '
                        + time
                        + '\n'
                        + 'Total number of samples: '
                        + str(len(LIST_ANNO))
                        + '\n'
                        + '# No Data: '
                        + str(len(list_nodata))
                        + '\n'
                        + '# Non-relevant: '
                        + str(len(list_irrelevant))
                        + '\n'
                        + '# Greater than 0.1% AF: '
                        + str(len(list_highfreq))
                        + '\n'
                        + '# Not expressed in brain: '
                        + str(len(list_notexpressedbrain))
                        + '\n'
                        + '# Candidates: '
                        + str(len(list_candidate))
                        + '\n'
                        + '# Candidate SNV:'
                        + str(len(list_snv))
                        + '\n'
                        + '# Candidate deletions:'
                        + str(len(list_del))
                        + '\n'
                        + '# Candidate insertions:'
                        + str(len(list_ins))
                        )
        summaryfile.close()
    elif option == 7:
        #Import annotations if not already in program:
        if len(LIST_CANDIDATE) == 0:
            filename = input('Please enter the name of the file being imported : ')
            LIST_CANDIDATE = va.importanno(filename)
        getsequences(LIST_CANDIDATE, BASEPATH)
    elif option == 8:
        check = input('Are you sure you want to exit? (y/n)')
        if check == 'y':
            EXIT_PROGRAM = True
        else:
            print('Returning to menu' + '\n')

print ("Terminating program...")
