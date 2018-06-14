from argparse import ArgumentParser
import re

#get the input and output files
parser = ArgumentParser()
parser.add_argument("--orfinput", dest="infile", help="input file with orf counts")
parser.add_argument("--gtf", dest="gtf", help="gtf file with mapping from orf to feature")
parser.add_argument("--feature",dest="feat", help="feature to enumerate from the gtf and orf counts")
parser.add_argument("--output", dest="outfile", help="output file with one annotation per row, each annotation carries counts from original feature")
args = parser.parse_args()

#open the files
infile = open(args.infile,"rU")
gtffile = open(args.gtf,"rU")
outfile = open(args.outfile,"w")
f = args.feat

#get the feature to orf mapping from the gtf
gtfdic = {}
for i in gtffile:
    row = i.strip("\n").split("\t")
    for j in row[-1].split(";"):
        if j != "":
            reg = re.search('^(\S+) \"(.*)\"$',j)
            if reg.group(1) == "gene_id":
                gtfdic[reg.group(2)]={}
                orf = reg.group(2)
            else:
                gtfdic[orf][reg.group(1)]=reg.group(2)

#dictionary to store unique features
featdic = {}
featlist = []

#parse input, storing counts as list and total length of each feature
for i in infile:
    if i[0] == "#" or re.match("Geneid",i):
        outfile.write(i)
    else:
        row = i.strip("\n").split("\t")
        geneid = row[0]
        contigs = row[1]
        starts = row[2]
        ends = row[3]
        strand = row[4]
        length = row[5]
        counts = row[6:]
        if f in gtfdic[geneid]:
            features = gtfdic[geneid][f].split(",")
            for j in features:
                if j in featdic:
                    featdic[j]["contigs"]+=";"+contigs
                    featdic[j]["starts"]+=";"+starts
                    featdic[j]["ends"]+=";"+ends
                    featdic[j]["strand"]+=";"+strand
                    featdic[j]["length"]+=int(length)
                    for k in range(len(featdic[j]["counts"])):
                        featdic[j]["counts"][k] += int(counts[k])
                else:
                    featdic[j]={}
                    featlist.append(j)
                    featdic[j]["contigs"]=contigs
                    featdic[j]["starts"]=starts
                    featdic[j]["ends"]=ends
                    featdic[j]["strand"]=strand
                    featdic[j]["length"]=int(length)
                    featdic[j]["counts"]=[int(x) for x in counts]
             
#write the parsed data to file
for i in featlist:
    outfile.write("\t".join([
        i,
        featdic[i]["contigs"],
        featdic[i]["starts"],
        featdic[i]["ends"],
        featdic[i]["strand"],
        str(featdic[i]["length"]),
        "\t".join([str(x) for x in featdic[i]["counts"]])
    ])+"\n")
outfile.close()
        
        
