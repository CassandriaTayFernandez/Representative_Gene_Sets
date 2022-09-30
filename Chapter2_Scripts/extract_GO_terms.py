#Input is the .tsv file produced by InterProScan
from glob import glob
from collections import defaultdict, Counter
annot_dict = defaultdict(list)
for x in glob('*TSV'):
    for line in open(x):
        ll = line.split('\t')
        # ['tripr.MilvusB.gnm2.ann1.mRNA27709_frame=1', '0768ee41eeda5f7c7526c231165147ff', '737', 'Pfam', 'PF00481', 'Protein phosphatase 2C', '50', '158', '1.6E-6', 'T', '05-04-2022', 'IPR001932', 'PPM-type phosphatase domain\n']
        pfam = '_'.join(ll[4:6])
        ll[0] = ll[0].replace('_frame=1','')
        annot_dict[ll[0]].append(pfam)

for line in open('Orthogroups.tsv'):
    ll = line.rstrip().split('\t')
    ortho = ll[0]
    count_dict = Counter()
    total = 0
    for spec in ll[1:]:
        genes = [x.strip() for x in spec.split(',')]
        for g in genes:
            total += 1
            if g in annot_dict:
                for pfam in annot_dict[g]:
                    count_dict[pfam] += 1
    group_string = []
    for x in count_dict:
        counts = count_dict[x]
        x = x.replace(' ','_')
        group_string.append('%s:%s'%(x, counts))
    print(ortho + '\t' + str(total) + '\t'+ '\t'.join(group_string))
