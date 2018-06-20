'''
classes and utility functions for pipeline_metaenumerate.py

'''

import os
import PipelineMetaAssemblyKit


'''
function to build bowtie database from input fasta
'''
def buildBowtieDB(infile,outfile,params):
    bcall = ["bowtie2-build {} {} -f --large-index".format(infile,outfile)]
    bcall.append("--threads {}".format(params["BowtieDB_threads"]))
    bcall.append("--seed {}".format(params["BowtieDB_seed"]))
    return(" ".join(bcall))


'''
function to call featureCount on bam files versus a GTF
'''
def countFeatures(feature,gtf,pairedness,output,inputs,params):
    cfcall=["featureCounts -a {} -o {} -g {} -t {}".format(gtf,output,feature,params["featureCounts_t"])]
    cfcall.append("-F {}".format(params["featureCounts_f"]))
    if params["featureCounts_o"] == "true":
        cfcall.append("-O")
    cfcall.append("--minOverlap {}".format(params["featureCounts_minoverlap"]))
    cfcall.append("--fracOverlap {}".format(params["featureCounts_fracoverlap"]))
    cfcall.append("--fracOverlapFeature {}".format(params["featureCounts_fracoverlapfeature"]))
    if params["featureCounts_largestoverlap"] == "true":
        cfcall.append("--largestOverlap")
    if params["featureCounts_readextension5"] != "false":
        cfcall.append("--readExtension5 {}".format(params["featureCounts_readextension5"]))
    if params["featureCounts_readextension3"] != "false":
        cfcall.append("--readExtension3 {}".format(params["featureCounts_readextension3"]))
    if params["featureCounts_read2pos"] != "false":
        cfcall.append("--read2pos {}".format(params["featureCounts_read2pos"]))
    if params["featureCounts_m"] == "true":
        cfcall.append("-M")
    if params["featureCounts_fraction"] == "true":
        cfcall.append("--fraction")
    cfcall.append("-s {}".format(params["featureCounts_s"]))
    if pairedness == True:
        if params["featureCounts_p"] == "true":
            cfcall.append("-P")
        if params["featureCounts_b"] == "true":
            cfcall.append("-B")
        if params["featureCounts_c"] == "true":
            cfcall.append("-C")
    cfcall.append("-T {}".format(params["featureCounts_threads"]))
    if params["featureCounts_v"] == "true":
        cfcall.append("--verbose")
    cfcall.append(" ".join(inputs))
    return(" ".join(cfcall))
