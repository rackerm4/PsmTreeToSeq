from setuptools import setup


setup(
    name='PsmTreeToSeq',
    description='Pipeline using DendroPy & Seq-Gen: Simulates pyhlogentic trees with protracted speciation model, '
                'using Seq-Gen to simulate the evolution of nucleotide or amino acid sequences along those '
                'phylogenies.',
    version='0.0.1',
    author = 'Chris P. Bacon',
    url = 'https://github.com/rackerm4/PsmTreeToSeq',
    long_description=open('README.md').read(),
    install_requires=[
            'dendropy>=4.4.0',
            'future==0.18.2',
            'PyYAML==5.3.1',
            'numpy==1.19.1'
      ],
)