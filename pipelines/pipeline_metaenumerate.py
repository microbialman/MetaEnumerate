"""
=============================
Metagenome enumeration pipeline
=============================

:Author: Matt Jackson
:Release: $Id$
:Date: |today|
:Tags: Python

Takes filtered reads from metagenomic data and enumerates annotations per sample from pre-defined ORFs within the data.

Overview
========
  
Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file. 


Input
-----

Reads
+++++

The pre-filtered reads to be aligned are imported by placing files are linking to files in the
:term:`working directory`.

The default file format assumes the following convention, with one file/paired file per sample to be enumerated:

   <filename>.<suffix>

The ``suffix`` determines the file type.  The
following suffixes/file types are possible:

fasta,fastq

fasta.1,fastq.1

fasta.gz,fastq.qz

fasta.1.gz,fastq.1.gz

Where unnumbered files relate to single end reads or inter-leaved paired-end files.
Interleaved files should be detected automatically.
Files containg 1 in the suffix will be assumed paired end and should be in the same 
directory as the read 2 files, which share the same file name except for the read number.

The location of the contigs and ORF annotation files must be provided in the pipeline.ini file 


.. note::

   Quality scores need to be of the same scale for all input files. Thus it might be
   difficult to mix different formats.

Optional inputs
+++++++++++++++

Requirements
------------

On top of the default CGAT setup and the MetaAssemblyKit pipeline, the pipeline requires the following
software to be in the path:

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Version*          |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
|bowtie2             | 2.3.0             | Mapping of reads to contigs                    |
+--------------------+-------------------+------------------------------------------------+

Pipeline output
===============

The main outputs are tables of counts of functional and taxonomic annotations per sample.

Glossary
========

.. glossary::

Code
====

"""
#load modules
from ruffus import *
import os, re
import sys
import subprocess


###################################################
###################################################
###################################################
# Pipeline configuration
###################################################
# load options from the config file
import CGATPipelines.Pipeline as P
P.getParameters(
       ["%s.ini" % __file__[:-len(".py")],
       "../pipeline.ini",
       "pipeline.ini" ] )
PARAMS = P.PARAMS

#add PipelineMetaAssemblyKit
import PipelineMetaAssemblyKit
import PipelineMetaFilter
import PipelineMetaEnumerate


#get all files within the directory to process
SEQUENCEFILES = ("*.fasta", "*.fasta.gz", "*.fasta.1.gz", "*.fasta.1",
                 "*.fna", "*.fna.gz", "*.fna.1.gz", "*.fna.1",
                 "*.fa", "*.fa.gz", "*.fa.1.gz", "*.fa.1", 
                 "*.fastq", "*.fastq.gz", "*.fastq.1.gz","*.fastq.1")

SEQUENCEFILES_REGEX = regex(
    r"(\S+).(fasta$|fasta.gz|fasta.1.gz|fasta.1|fna$|fna.gz|fna.1.gz|fna.1|fa$|fa.gz|fa.1.gz|fa.1|fastq$|fastq.gz|fastq.1.gz|fastq.1)")


#######################################################
# Make a bowtie indexed database from the contigs
#######################################################
@follows(mkdir("contig_database.dir"))
@transform(PARAMS["General_contig_file"],regex(r"(\S+)/(\S+).(fasta|fa|fna)($|.gz$)"),r"contig_database.dir/\2.bowtie.1.bt2")
def makeBowtieDb(infile,outfile):
    #check input is a signle end fasta as appropriate
    seqdat = PipelineMetaAssemblyKit.SequencingData(infile)
    if seqdat.paired == False and seqdat.fileformat == "fasta":
        #call to bowtie2-build
        job_memory = str(PARAMS["BowtieDB_memory"])+"G"
        job_threads = PARAMS["BowtieDB_threads"]
        statement = PipelineMetaEnumerate.buildBowtieDB(infile,outfile.replace(".1.bt2",""),PARAMS)
        P.run()
    else:
        print("Bowtie database can only be constructed from a single end FASTA file. Please check contigfile in parameters.")


#######################################################
# Map each samples reads against the contig database
#######################################################
@follows(makeBowtieDb)
@follows(mkdir("sample_mappings.dir"))
@transform(SEQUENCEFILES,SEQUENCEFILES_REGEX,r"sample_mappings.dir/\1/\1.mapped.bam")
def mapSample(infile,outfile):
    #set the mapping db to the newly created index
    PARAMS["Bowtie_genome_db"]="contig_database.dir/"+re.search(r"(\S+)/(\S+).(fasta|fa|fna)($|.gz$)",PARAMS["General_contig_file"]).group(2)+".bowtie"
    #use bowtie call from the MetaFilter pipeline
    job_threads = PARAMS["Bowtie_threads"]
    job_memory = str(PARAMS["Bowtie_memory"])+"G"
    seqdat = PipelineMetaAssemblyKit.SequencingData(infile)
    bowtie = PipelineMetaFilter.Bowtie2(seqdat,outfile,PARAMS)
    #need to reset the working directory in the bowtie function as it is running on files in one directory (opposed to multiple in MetaFilter pipeline)
    bowtie.indir = ""
    statementlist = []
    #remove all comments from read names in files (trimming can add comments making non-matching readnames in pairs)
    statementlist.append(bowtie.cleanNames())
    #directory for output
    statementlist.append("mkdir -p {}".format(os.path.dirname(outfile)))
    #call to bowtie
    statementlist.append(bowtie.build())
    #convert sam to bam
    statementlist.append("samtools view -bS {} > {}".format(outfile.replace(".bam",".sam"),outfile))
    #remove the sam file
    statementlist.append("rm {}".format(outfile.replace(".bam",".sam")))
    statement = " && ".join(statementlist)
    P.run()


#################################################################
# Count reads mapping to each ORF annotation using featureCount
#################################################################
@follows(mapSample)
@follows(mkdir("feature_counts.dir"))
@merge(mapSample,"feature_counts.dir/orf_counts.txt")
def countFeatures(infiles,outfile):
    #get pairdness from first sample
    paired = True
    num_paired = subprocess.check_output(["samtools","view","-c","-f 1","{}".format(os.getcwd()+"/"+infiles[0])]).decode(sys.stdout.encoding)
    if int(num_paired.strip("\n")) == 0:
        paired = False
    #generate counts per orf across all samples
    job_threads = PARAMS["featureCounts_threads"]
    statementlist = [PipelineMetaEnumerate.countFeatures("gene_id",PARAMS["General_gtf_file"],paired,outfile,infiles,PARAMS)]
    #generate counts for other features specified in params
    for i in PARAMS["General_feature_list"].split(","):
        statementlist.append(PipelineMetaEnumerate.countFeatures(i,PARAMS["General_gtf_file"],paired,outfile.replace("orf",i),infiles,PARAMS))
    statement = " && ".join(statementlist)
    P.run()
    
@follows(countFeatures)
def full():
    pass


if __name__ == "__main__":
    if sys.argv[1] == "plot":
        pipeline_printout_graph("test.pdf", "pdf", [full], no_key_legend=True,
                                size=(4, 4),
                                user_colour_scheme = {"colour_scheme_index": 1})
    else:
        sys.exit(P.main(sys.argv))
