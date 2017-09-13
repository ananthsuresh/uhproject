import sys
import numpy as np
import pandas as pd
import scipy.stats as stats
import pylab as pl
import matplotlib.pyplot as plt
from fuzzysearch import find_near_matches
import editdistance
import subprocess
import os
# Global Variables
filename = sys.argv[1]
directory = "/users/ananth/desktop/samples1/"
longestClippedRead = ""
totalReadCount = 0
singleFlankMutant = 0
bothFlankMutant = 0
singleFlankNonMutant = 0
bothFlankNonMutant = 0
noFlank = 0

def extractReadsFastQ(filename):
    with open(filename , 'r') as f:
        f.readline()
        linenumber = 0
        for line in f:
            if (linenumber % 4 == 0):
                yield line
            linenumber += 1
    f.close()

def QC(minQuality,minLength):
    global filename
    print("Filtering reads by Quality")
    os.system("sickle se -f " + filename + ".fastq -t sanger -o " + filename + "_qc.fastq -q " + str(minQuality) + " -l " + str(minLength))
    filename = sys.argv[1] + "_qc"

def extractFLT3():
    global filename
    print("Aligning to HG19 - FLT3")
    os.system("bowtie2 -x " + directory + "hg19_minusFLT3 " + filename + ".fastq -S " + filename + "_HG19-FLT3.sam" )
    print("Converting to BAM")
    os.system("samtools view -bS " + filename + "_HG19-FLT3.sam > " + filename + "_HG19-FLT3.bam")
    print("Extracting badly mapped reads(FLT3)")
    os.system("samtools view -h -b -q 30 -U " + filename + "_onlyFLT3.bam " + filename + "_HG19-FLT3.bam -o aboveq.bam")
    print("Sorting BAM file")
    os.system("samtools sort " + filename + "_onlyFLT3.bam -o " + filename + "_onlyFLT3_sorted.bam")
    print("Converting to FASTQ")
    os.system("samtools bam2fq " + filename + "_onlyFLT3_sorted.bam > " + filename + "_onlyFLT3.fastq")
    filename += "_onlyFLT3"

def getConsensusSequence():
    global filename
    makeFasta(longestClippedRead)
    print("Indexing Longest Read")
    os.system("bowtie2-build " + filename + "_longestRead.fa " + filename + "_longestRead")
    print("Aligning all mutated reads to longest mutated read")
    os.system("bowtie2 --local -x " + filename + "_longestRead " +  filename + ".fastq -S " + filename + "_alnToLongestRead.sam")
    print("Converting to BAM")
    os.system("samtools view -bS " + filename + "_alnToLongestRead.sam > " + filename + "_alnToLongestRead.bam")
    print("Sorting BAM file")
    os.system("samtools sort " + filename + "_alnToLongestRead.bam -o " + filename + "_alnToLongestRead_sorted.bam")
    print("Indexing BAM File")
    os.system("samtools index " + filename + "_alnToLongestRead_sorted.bam " + filename + "_alnToLongestRead_sorted.bai")
    print("Getting Consensus Sequence")
    os.system("samtools mpileup -uf " + filename + "_longestRead.fa " + filename + "_alnToLongestRead_sorted.bam | bcftools call -c | vcfutils.pl vcf2fq > " + filename + "_consensus.fastq")

def slidingWindow(sequence,winSize,stepSize=1):
    # Pre-compute number of chunks to emit
    numOfChunks = ((len(sequence)-winSize) // stepSize)+1

    # Do the work
    for i in range(0,numOfChunks*stepSize,stepSize):
        yield sequence[i:i+winSize]

def makeFastQ(mutatedReads):
    global filename
    index = 0
    mutations = open(filename + "_mutated.fastq", 'w')
    for read in mutatedReads:
        mutations.write("@" + str(index) + "\n")
        read = read.rstrip()
        mutations.write(read + "\n")
        mutations.write("+\n")
        quality = "~" * len(read)
        mutations.write(quality + "\n")
        index += 1
    mutations.close()
    filename += "_mutated"
def makeFasta(longestClippedRead):
    longestRead = open(filename + "_longestRead.fa", 'w')
    longestRead.write(">" + sys.argv[1] + "\n")
    longestRead.write(longestClippedRead)
    longestRead.close()

def plotHistogram(mutatedLengths, unmutatedLengths):
    unmutated = sorted(unmutatedLengths)
    mutated = sorted(mutatedLengths)
    pl.hist(unmutated,normed=False,range = [100,200], label = 'unmutated')
    pl.hist(mutated,normed=False, range = [100,200],  label = 'mutated')
    pl.legend(loc='upper right')
    pl.title("Histogram of Length -" + sys.argv[1])
    pl.savefig(sys.argv[1] + 'length.png')
    pl.close()

def plotStats():
    raw_data = {'Class': ['Single Flank', 'Both Flank'],
            'Mutant': [singleFlankMutant, bothFlankMutant],
            'Non Mutant': [singleFlankNonMutant, bothFlankNonMutant]}
    df = pd.DataFrame(raw_data, columns = ['Class', 'Mutant', 'Non Mutant'])
        # Create the general blog and the "subplots" i.e. the bars
    f, ax1 = plt.subplots(1, figsize=(10,5))

    # Set the bar width
    bar_width = 0.75

    # positions of the left bar-boundaries
    bar_l = [i for i in range(len(df['Mutant']))]

    # positions of the x-axis ticks (center of the bars as bar labels)
    tick_pos = [i for i in bar_l]

    # Create a bar plot, in position bar_1
    ax1.bar(bar_l,
            # using the pre_score data
            df['Mutant'],
            # set the width
            width=bar_width,
            # with the label pre score
            label='Mutant',
            # with alpha 0.5
            alpha=0.5)

    # Create a bar plot, in position bar_1
    ax1.bar(bar_l,
            # using the mid_score data
            df['Non Mutant'],
            # set the width
            width=bar_width,
            # with pre_score on the bottom
            bottom=df['Mutant'],
            # with the label mid score
            label='Non Mutant',
            # with alpha 0.5
            alpha=0.5)


    # set the x ticks with names
    plt.xticks(tick_pos, df['Class'])

    # Set the label and legends
    ax1.set_ylabel("Frequency")
    ax1.set_xlabel("Type")
    plt.legend(loc='upper right')

    # Set a buffer around the edge
    plt.xlim([min(tick_pos)-bar_width, max(tick_pos)+bar_width])
    plt.savefig(sys.argv[1] + 'stats.png')
    plt.close()

def detectMutation():
    print("Searching for mutant reads...")
    global longestClippedRead
    global totalReadCount
    global singleFlankMutant
    global bothFlankNonMutant
    global singleFlankNonMutant
    global bothFlankMutant
    global noFlank

    maxMismatch = 1
    rightFlank = 'TTCATACCTA'
    leftFlank = 'CTTACCAAAC'
    insideReference = 'TCTAAATTTTCTCTTGGAAACTCCCATTTGAGATCATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGGTCACCTGTACCATCTGTAGCTGGCT'
    lengthThreshold = len(insideReference) + 25
    mutatedReads = []
    unmutatedReads = []
    mutatedLengths = []
    unmutatedLengths = []
    maxLength = 0
    stepSize = 1
    windowSize = 5
    reads = extractReadsFastQ(filename + ".fastq")

    for read in reads:
        totalReadCount += 1
        # print(totalReadCount)
        #looks for match on left and right flank with mismatch
        leftMatches = find_near_matches(leftFlank,read, max_l_dist=maxMismatch)
        rightMatches = find_near_matches(rightFlank,read,max_l_dist=maxMismatch)
        #if there are no matches on either side, go on to next read
        if not leftMatches and not rightMatches:
            noFlank += 1
            continue
        #if only left matches, cut the left out and calculate edit distance of what's left
        if leftMatches and not rightMatches:
            #cuts out left flanking region
            clippedRead = read[leftMatches[0].start + len(leftFlank):]
            if(len(clippedRead) > maxLength):
                maxLength = len(clippedRead)
                longestClippedRead = clippedRead

            if(len(clippedRead) > len(insideReference)):
                mutatedReads.append(clippedRead)
                singleFlankMutant += 1
                continue
            chunks = slidingWindow(clippedRead, windowSize, stepSize)
            i = 0
            mutated = False
            for chunk in chunks:
                score = editdistance.eval(chunk,insideReference[i * stepSize:len(chunk) + (i * stepSize)])
                if(score >= len(chunk)):
                    mutatedReads.append(clippedRead)
                    singleFlankMutant += 1
                    mutated = True
                    break #need to continue here
                i += 1
            if not mutated:
                unmutatedReads.append(clippedRead)
                singleFlankNonMutant += 1
                continue
            continue
        if rightMatches and not leftMatches:
            clippedRead = read[:rightMatches[0].start]
            if(len(clippedRead) > maxLength):
                maxLength = len(clippedRead)
                longestClippedRead = clippedRead
            if(len(clippedRead) > len(insideReference)):
                mutatedReads.append(clippedRead)
                singleFlankMutant += 1
                continue
            chunks = slidingWindow(clippedRead, windowSize, stepSize)
            i = 0
            mutated = False
            for chunk in chunks:
                if(i == 0):
                    score = editdistance.eval(chunk,insideReference[-(len(chunk) + (i * stepSize)):])
                else:
                    score = editdistance.eval(chunk,insideReference[-(len(chunk) + (i * stepSize)):-(i * stepSize)])
                if(score >= len(chunk)):
                    mutatedReads.append(clippedRead)
                    singleFlankMutant += 1
                    mutated = True
                    break #need to continue here
                i += 1
            if not mutated:
                unmutatedReads.append(clippedRead)
                singleFlankNonMutant += 1
                continue
            continue
        if rightMatches and leftMatches:
            clippedRead = read[leftMatches[0].start + len(leftFlank): rightMatches[0].start]
            if(len(clippedRead) > maxLength):
                maxLength = len(clippedRead)
                longestClippedRead = clippedRead
            if(len(clippedRead) > lengthThreshold):
                mutatedReads.append(clippedRead)
                mutatedLengths.append(len(clippedRead))
                bothFlankMutant += 1
            else:
                unmutatedReads.append(clippedRead)
                unmutatedLengths.append(len(clippedRead))
                bothFlankNonMutant += 1
                continue

    makeFastQ(mutatedReads)
    plotHistogram(mutatedLengths, unmutatedLengths)
    plotStats()

if __name__ == '__main__':
    os.chdir(filename)
    QC(20,80)
    extractFLT3()
    detectMutation()
    getConsensusSequence()
