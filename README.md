# flyingoctopus
Variantannotation

This program is a wrapper for myvariant.info and Ensembl's VEP REST API.

Allows for annotation of mutations straight from .vcf files.

## Getting Started

1. Install dependencies

    With virtualenv run `make venv-install`

    Without run `make install`

2. Run

    With virtualenv run
    ```
    source venv/bin/activate
    python VariantAnnotation.py
    ```

    Without run `python VariantAnnotation.py`