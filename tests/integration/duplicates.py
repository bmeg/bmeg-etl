import subprocess


def test_g2p():
    """ check all gids (except outputs/g2p/Deadletter.Vertex.json )"""
    cmd = 'cat outputs/g2p/Allele.Vertex.json outputs/g2p/Compound.Vertex.json outputs/g2p/G2PAssociation.Vertex.json outputs/g2p/MinimalAllele.Vertex.json outputs/g2p/Phenotype.Vertex.json | jq .gid | sort | uniq -d'
    output = subprocess.check_output(cmd, shell=True).decode()
    assert len(output.split()) == 0, 'G2P transformer should not emit duplicate gid\n{}'.format(output)


def test_compound():
    """ check all gids"""
    cmd = 'cat outputs/compound/normalized.Compound.Vertex.json | jq .gid | sort | uniq -d'
    output = subprocess.check_output(cmd, shell=True).decode()
    assert len(output.split()) == 0, 'Compound transformer should not emit duplicate gid\n{}'.format(output)


def test_phenotype():
    """ check all gids"""
    cmd = 'cat outputs/phenotype/normalized.Phenotype.Vertex.json | jq .gid | sort | uniq -d'
    output = subprocess.check_output(cmd, shell=True).decode()
    assert len(output.split()) == 0, 'Phenotype transformer should not emit duplicate gid\n{}'.format(output)
