# flyingoctopus
Variantannotation

This program is a wrapper for myvariant.info and Ensembl's VEP REST API.

Allows for annotation of mutations straight from .vcf files.

splitfile branch - purpose of this branch is to split the va.py file into multiple files, organized by the annotation source. 

Essentially, the import and export functions will be moved to Variantannotation, while the annotation and parsing function will be separated into their respective files.
