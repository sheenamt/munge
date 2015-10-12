==========================================================
Scripts for munging output of NGS pipeline
==========================================================

This project provides a program with a command-line interface for
parsing Next Generation Sequencing data.

.. contents:: Table of Contents

dependencies
============

* Python 2.7.x
* A UNIX-like operating system (Linux, OS X). Not tested on Windows.

installation
============

Clone the project from the git repository::

    cd ~/src
    git clone git@bitbucket.org:uwlabmed/munge.git
    cd munge

Now installation can be performed using the install script provided.
This will default to install at /home/genetics unless an install path is provided.

    sudo ./install_munge /install/path

This script does a clean install.


execution
=========

The ``munge`` script provides the different scripts used to process
data output from the pipeline. Note that for development, it is convenient
to run ``munge`` from within the project directory by specifying the
relative path to the script::

    % ./munge

Commands are constructed as follows. Every command starts with the
name of the script, followed by an "action" followed by a series of
required or optional "arguments". The name of the script, the action,
and options and their arguments are entered on the command line
separated by spaces. Help text is available for both the ``munge``
script and individual actions using the ``-h`` or ``--help`` options::

    % munge -h
    usage: munge [-h] [-V] [-v] [-q]

             {help,xlsmaker,rename_hiseq,sample_crawler,}...
Utilities for the munge scripts

positional arguments:
  {help,xlsmaker,rename_hiseq,control_parser,variant_crawler,
  freq_creator,rename_miseq,db_annotation,quality_metrics,
  getpfx,combined_cnv,combined_output,annovar_bed_parser,
  qc_variants,combined_pindel,summary}

    help                Detailed help for actions using `help <action>`
    xlsmaker            Create xls workbook from all output files
    rename_hiseq        Rename and compress HiSeq files.
    control_parser      Compare quality control variants to OPX-240 output to
                        check quality of run
    variant_crawler     Create annovar file from Clinical variants csv
    freq_creator        Calculate tallies of variants and write anovar output
    rename_miseq        Rename MiSeq files for pipeline processing
    db_annotation       Create annotation of all variants in db (or only from
                        GATK)
    quality_metrics     Parse picard and CNV output to create quality metrics
                        file
    getpfx              Get prefixes files (PFX.[12].fastq.gz) for running
                        pipeline.
    combined_cnv        Crawl analysis files to create one analysis file with
                        all info
    combined_output     Crawl analysis files to create one analysis file with
                        all info
    annovar_bed_parser  Filter a file of genomic positions given ranges of
                        start positions
    qc_variants         Parse variant files from pipeline, 1000G, and Complete
                        Genomics to create QC Variant file
    combined_pindel     Crawl analysis files to create one analysis file with
                        all info
    summary             Summarize output from Annovar and EVS

optional arguments:
  -h, --help            show this help message and exit
  -V, --version         Print the version number and exit
  -v, --verbose         Increase verbosity of screen output (eg, -v is
                        verbose, -vv more so)
  -q, --quiet           Suppress output

Help text for an individual action is available by including the name
of the action::
    % munge getpfx -h
    usage: munge getpfx [-h] [-s SEPARATOR] datadir

    Get prefixes files (PFX.[12].fastq.gz) for running pipeline.

    Usage:
        munge getpfx /path/to/data

    positional arguments:
        datadir               Path to directory containing fastq files.

    optional arguments:
        -h, --help            show this help message and exit
        -s SEPARATOR, --separator SEPARATOR
                        separator for list of prefixes

versions
========

We use abbrevited git sha hashes to identify the software version::

    % ./munge -V
    0309.004ecac

unit tests
==========

Unit tests are implemented using the ``unittest`` module in the Python
standard library. The ``tests`` subdirectory is itself a Python
package that imports the local version (ie, the version in the project
directory, not the version installed to the system) of the ``munge``
package. All unit tests can be run like this::

    munge % ./testall
    ........................
    ----------------------------------------------------------------------
    Ran 24 tests in 0.155s

    OK

A single unit test can be run by referring to a specific module,
class, or method within the ``tests`` package using dot notation::

    munge % ./testone tests.test_subcommands.TestQCVariants
    .
    ----------------------------------------------------------------------
    Ran 1 test in 0.004s

    OK

