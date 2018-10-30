

# BMEG Neo4j

## Scope

* Evaluate Neo4j support for BMEG's graph oriented `vertex` and `edge` types.
* Evaluate Neo4j support for BMEG's document model, specifically json queries

## Out of Scope

* "Bake-off" with mongo
* Translation and strategies for "grip to cypher" i.e grip queries query translation (all queries were hand crafted)


## neo4j setup

```
alias cs='cypher-shell'
export TERM=xterm
```

/etc/neo4j/neo4j.conf

```

# *****************
# local changes
# ****************
dbms.security.procedures.unrestricted=apoc.*
apoc.export.file.enabled=true
apoc.import.file.enabled=true
apoc.import.file.use_neo4j_config=true
dbms.security.allow_csv_import_from_file_urls=true
```


cd /var/lib/neo4j/plugins/
```
sudo wget https://github.com/neo4j-contrib/neo4j-apoc-procedures/releases/download/3.4.0.3/apoc-3.4.0.3-all.jar
sudo chown neo4j apoc-3.4.0.3-all.jar
sudo chgrp adm apoc-3.4.0.3-all.jar
```


Use cases:
* add papers associated with genes
