from setuptools import setup, find_packages, Extension
from codecs import open
from os import path


triemodule = Extension('trie',
                       sources = ['trie/triemodule.c', 'trie/trie.c'])

# Get the long description from the README file
here = path.abspath(path.dirname(__file__))
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    # TODO: specify dependencies
    name = 'MultiGuidescan',
    version = '1.0',
    author= 'Alexendar Perez & Yuri Pritykin',
    author_email= 'pereza1@mskcc.org & pritykiy@mskcc.org',
    description = 'Tools to create and interface genome-wide CRISPR'
                  ' guideRNA databases',
    long_description = long_description,

    install_requires = [
        'biopython>=1.66','pysam==0.8.3','pyfaidx==0.4.7.1','bx-python==0.7.3'
    ],
    dependency_links = ['https://github.com/bxlab/bx-python.git'
                        ],
    packages = find_packages(),
    ext_modules = [triemodule],
    package_data = {'guidescan' : ['annotation_bed/*/*.bed','Rule_Set_2_scoring/saved_models/V3_model_*.pickle','CFD_scoring/*.pkl']},
    #test_suite = 'tests',
    entry_points={
        'console_scripts': [
            'guidescan_processer = guidescan.processer:main',
            'guidescan_bamdata = guidescan.bamdata:main',
            'guidescan_guidequery = guidescan.guidequery:main',
            'guidescan_cutting_efficiency_processer = guidescan.cutting_efficiency_processer:main',
            'guidescan_cutting_specificity_processer = guidescan.cutting_specificity_processer:main'
        ],
    }
)
