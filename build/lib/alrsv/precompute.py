import json

from Bio import SeqIO


def process_reference():
    cov = list(SeqIO.parse("sequence.gb", "genbank"))[0]

    genes = {}

    for f in cov.features:
        if f.type == 'gene':
            gene = f.qualifiers['gene'][0]
            start = int(f.location.start)
            end = int(f.location.end)
            genes[gene] = [start, end]

    with open('sars_cov_2.py', 'w') as f:
        f.write('genes = {}\n\nseq = \'{}\''.format(genes, cov.seq))


#def get_amplicons():
#    cov = list(SeqIO.parse("sequence.gb", "genbank"))[0]
#    with open('nCoV-2019.insert.bed', 'r') as f:
#        inserts = [l.split('\t') for l in f.read().split('\n')[:-1]]
#
#    for i in range(len(inserts)):
#        insert = inserts[i]
#        section = cov.seq[int(insert[1]):int(insert[2])]
#        gc = (section.count('G') + section.count('C')) / len(section)
#        inserts[i].append(gc)
#
#    with open('artic_amplicons.py', 'w') as f:
#        f.write('inserts = {}'.format(inserts))

def get_clades():
    import json
    from os import getenv, listdir, mkdir, system
    from os.path import isfile, join

    def parse_mut(mut):
        mut = mut.upper()
        return mut

    data_paths = [
        './data/clades/',
    ]

    clades = {}
    mutations = {}
    for data_path in data_paths:
        fns = [f for f in listdir(data_path) if isfile(join(data_path, f)) and f.endswith('.json')]
        fns.sort()
        for fn in fns:
            with open('{}/{}'.format(data_path, fn), 'r') as f:
                voc = json.loads(f.read())
            label = voc['label']
            clades[label] = [parse_mut(mut) for mut in voc['sites']]
    with open('clades.py', 'w') as f:
         f.write('clades = {}'.format(clades))
    vocs = list(clades.keys())
    
    for voc in vocs:
        for mut_name in clades[voc]:
            if mut_name not in mutations:
                mutations[mut_name] = {voc: 0 for voc in vocs}
            mutations[mut_name][voc] = 1
    with open('mutations.py', 'w') as f:
         f.write('mutations = {}'.format(mutations))


if __name__ == '__main__':
    # process_reference()
    # get_amplicons()
    get_clades()
