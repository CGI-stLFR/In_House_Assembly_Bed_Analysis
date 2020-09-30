from subprocess import Popen, PIPE
import sys

bed_file = sys.argv[1]
vcf = sys.argv[2]

with open(bed_file, "r") as bed:
    for line in bed:
        parts = line.strip().split()
        
        reg = f"{parts[0]}:{parts[1]}-{parts[2]}"
        
        vcf_entries = Popen(["tabix", vcf, reg], stdout=PIPE, stderr=PIPE)
        vcfout, vcferr = vcf_entries.communicate()
        entries = vcfout.decode('utf-8').split('\n')
        variants = 0
        #variants_list = []
        for entry in entries:
            if entry.startswith("#") or not entry:
                continue
            else:
                entry = entry.split("\t")
                if entry[11].split(":")[0] != entry[12].split(":")[0]:
                    variants += 1
                    #variants_list.append(entry)
        if variants:
            parts.extend([int(parts[2]) - int(parts[1]), variants, variants/((int(parts[2]) - int(parts[1])) / 1000)])
        else:
            parts.extend([int(parts[2]) - int(parts[1]), 0, 0])
        
        parts = [str(x) for x in parts]
        print('\t'.join(parts))
