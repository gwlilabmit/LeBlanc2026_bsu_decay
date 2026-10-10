### removes linker sequence from Illumina sequencing results ###



import numpy as np

import os, subprocess

import sys



def stepAlign(fqIN, refs, align_options):

    steps = len(refs)

    fq = fqIN

    for i in range(steps):

    

        alignedFile = fq[:fq.find('.')] + '_' + refs[i][:refs[i].find('/')] + '_' + align_options + '.align'

        unalignedFile = fq[:fq.find('.')] + '_' + refs[i][:refs[i].find('/')] + '_' + align_options + '.unalign'

        

        if align_options == 'v0m1':

            output=subprocess.Popen('bowtie -v 0 -m 1 --quiet /home/mtien/genomes/'+refs[i]+' '+fq+' '+alignedFile+' --un '+unalignedFile, shell=True, stdout=subprocess.PIPE)

            

        elif align_options == 'v1m1':

            output=subprocess.Popen('bowtie -v 1 -m 1 --quiet /home/mtien/genomes/'+refs[i]+' '+fq+' '+alignedFile+' --un '+unalignedFile, shell=True, stdout=subprocess.PIPE)

            

        elif align_options == 'v1m2k2':

            output=subprocess.Popen('bowtie -v 1 -m 2 -k 2 --quiet /home/mtien/genomes/'+refs[i]+' '+fq+' '+alignedFile+' --un '+unalignedFile, shell=True, stdout=subprocess.PIPE)

            

        elif align_options == 'v1k1':

            output=subprocess.Popen('bowtie -v 1 -k 1 /home/mtien/genomes/'+refs[i]+' '+fq+' '+alignedFile+' --un '+unalignedFile, shell=True, stdout=subprocess.PIPE)
            
        
        elif align_options == 'v2k1':

            output=subprocess.Popen('bowtie -v 2 -k 1 /home/mtien/genomes/'+refs[i]+' '+fq+' '+alignedFile+' --un '+unalignedFile, shell=True, stdout=subprocess.PIPE)
        

        output.wait()

	

        



if __name__ == '__main__':



    fqIN = sys.argv[1]

    sample_type = sys.argv[2]

    align_options = sys.argv[3]

    # v0m1, v1m1, v1m2k2, v1k1

    

    if sample_type=='bsub':

        stepAlign(fqIN, ['bsubtilis/bsub_chr'], align_options)

    elif sample_type=='ivt_spike':

        stepAlign(fqIN, ['ivt_spikein/IVT_spike/ivt_spike'], align_options)

    elif sample_type=='ivt4_spike':

        stepAlign(fqIN, ['ivt_spikein/ivt_ind_spike/ivt4_spike'], align_options)