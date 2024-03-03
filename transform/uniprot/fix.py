

def force_list(x):
    if isinstance(x, list):
        return x
    return [x]

def ensemblInfo(row):
    ensembl_gene = []
    ensembl_transcript = []
    ensembl_pseq = []
    for i in force_list(row['dbReference']):
        if i['-type'] == 'Ensembl':
            ensembl_transcript.append(i["-id"].split(".")[0])
            for p in force_list(i["property"]):
                if p["-value"].startswith("ENSG"):
                    ensembl_gene.append(p["-value"].split(".")[0])
                elif p['-value'].startswith('ENSP'):
                    ensembl_pseq.append( p["-value"].split(".")[0])
    o = []
    if 'feature' in row and len(ensembl_transcript) > 0:
        for f in force_list(row['feature']):
            c = {'transcript': ensembl_transcript[0]} #, 'sequence_id': ensembl_pseq} 
            description = f['-type']
            if '-description' in f:
                c[description.replace('-', '_').replace(' ', '_')] = f['-description']
            else:
                c['structure'] = description
            if 'location' in f:
                if 'begin' in f['location'] and 'end' in f['location']:
                    if '-position' in f['location']['begin']:
                        c['start'] = int(f['location']['begin']['-position'])
                    else:
                        c['start'] = None
                    if '-position' in f['location']['end']:
                        c['end'] = int(f['location']['end']['-position'])
                    else:
                        c['end'] = None
                elif 'position' in f['location'] and '-position' in f['location']['position']:
                    c['start'] = int(f['location']['position']['-position'])
                    c['end'] = int(f['location']['position']['-position'])
                else:
                    continue
                o.append(c)
    return o