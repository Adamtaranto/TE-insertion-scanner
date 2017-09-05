from setuptools import setup

pypi_classifiers = [
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3',
    "Development Status :: 4 - Beta",
    "Environment :: Console",
    "Operating System :: OS Independent",
    'Intended Audience :: Science/Research',
    'Natural Language :: English',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    "Topic :: Software Development :: Libraries :: Python Modules",
    'License :: OSI Approved :: MIT License',
]

install_requires = [
    'biopython>=1.70',
]

desc = """Find alignment signatures characteristic of transposon insertion sites."""

setup(name='tinscan',
      version='0.1.0',
      description=desc,
      url='https://github.com/Adamtaranto/TE-insertion-scanner',
      author='Adam Taranto',
      author_email='adam.taranto@anu.edu.au',
      license='MIT',
      packages=['tinscan'],
      classifiers=pypi_classifiers,
      keywords=["Transposon","TIR","MITE","TE","WGA","LASTZ","Whole genome alignment","insertion","transposition"],
      install_requires=install_requires,
      include_package_data=True,
      package_data={
        'src':[
             'src/LASTZ_genome_align.sh'],},
      zip_safe=False,
      entry_points={
        'console_scripts': [
            'tinscan-prep=tinscan.run_prep:main',
            'tinscan-align=tinscan.run_align:main',
            'tinscan-find=tinscan.run_scan:main',
        ],
    },
    )
