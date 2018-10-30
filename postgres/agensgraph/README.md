# AgensGraph

## Install it:
  https://hub.docker.com/r/bitnine/agensgraph/

## Load it:
  See [loader.py](loader.py) and [config.yml](config.yml)

##  Issues:

  * volume mapping of data directory to host not working
  https://github.com/bitnine-oss/agensgraph-docker#usage-docker
  " not a database cluster"
    * resolution - don't start that way - e.g. do not map to host volume /home/agens/AgensGraph/data
    ```
    volumes:
      - ./:/src
      # - ./postgres/agensgraph/data:/home/agens/AgensGraph/data    
    ```

  * by default the postgres instance in the supplied docker file `does not` connect to external hostss
    * resolution - log onto docker image, edit
       vi /home/agens/AgensGraph/data/pg_hba.conf
       vi /home/agens/AgensGraph/data/postgresql.conf
       ```
       pg_hba.conf
       host all all 0.0.0.0/0 trused
       postgresql.conf
       listen_addresses='*'
       ```

  * cypher weirdness
    * apparently the 52K field object is treated like a hash

  ```
  * given this structure
    MATCH (n:Expression {gid: 'Expression:gtex:GTEX-1117F-0226-SM-5GZZ7'} ) return jsonb_pretty(n) ;
                         jsonb_pretty
   --------------------------------------------------------
    {                                                     +
        "id": "GTEX-1117F-0226-SM-5GZZ7",                 +
        "gid": "Expression:gtex:GTEX-1117F-0226-SM-5GZZ7",+
        "method": "Illumina HiSeq",                       +
        "metric": "GENE_TPM",                             +
        "source": "gtex",                                 +
        "values": {                                       +
            "ENSG00000000003": 46.040000915527344,        +
            "ENSG00000000005": 16.579999923706055,        +
            "ENSG00000000419": 62.72999954223633,         +
            "ENSG00000000457": 7.584000110626221,         +


  # query works and row count ok, but no data?

  agens=# MATCH (n:Expression {gid: 'Expression:gtex:GTEX-1117F-0226-SM-5GZZ7'} ) return n.values.ENSG00000000003 ;
   ensg00000000003
  -----------------

  (1 row)

  # query works and row count ok, and data returned as expected.

  Time: 3.768 ms
  agens=# MATCH (n:Expression {gid: 'Expression:gtex:GTEX-1117F-0226-SM-5GZZ7'} ) return n.values['ENSG00000000003'] ;
         values
  --------------------
   46.040000915527344
  (1 row)

  # other properties are accessed with dot notation, as expected.

  MATCH (n:Expression {gid: 'Expression:gtex:GTEX-1117F-0226-SM-5GZZ7'} ) return n.metric ;
   metric
------------
 "GENE_TPM"
(1 row)

  ```


## json query performance

* without index

```

explain MATCH (n:Expression  ) where n.values['ENSG00000000003'] > 50  return count(*) ;
                                      QUERY PLAN
---------------------------------------------------------------------------------------
 Aggregate  (cost=241.84..241.85 rows=1 width=8)
   ->  Seq Scan on expression n  (cost=0.00..232.10 rows=3896 width=0)
         Filter: (properties.'values'::text['"ENSG00000000003"'::jsonb] > '50'::jsonb)


 agens=# MATCH (n:Expression  ) where n.values['ENSG00000000003'] > 50  return count(*) ;
  count
 -------
   1357
 (1 row)

 Time: 69202.875 ms
```

* create index on property

```

 CREATE PROPERTY INDEX expression_ENSG00000000003 ON Expression (values['ENSG00000000003']) ;
 Time: 45081.250 ms
```

* not shown in standard explain ?

```
 agens=# explain MATCH (n:Expression  ) where n.values['ENSG00000000003'] > 50  return count(*) ;
                                       QUERY PLAN
 ---------------------------------------------------------------------------------------
  Aggregate  (cost=241.84..241.85 rows=1 width=8)
    ->  Seq Scan on expression n  (cost=0.00..232.10 rows=3896 width=0)
          Filter: (properties.'values'::text['"ENSG00000000003"'::jsonb] > '50'::jsonb)
 (3 rows)

# but is shown in display graph vertex

# \dGv Expression;
              List of labels
   Graph    |    Name    |  Type  | Owner
------------+------------+--------+-------
 bmeg_graph | expression | vertex | agens
(1 row)

Vertex label "bmeg_graph.expression"
--
Property Indexes:
    "expression_ensg00000000003" btree (values[ENSG00000000003])
    "expression_gid_idx" UNIQUE btree (gid)
Inherits: bmeg_graph.ag_vertex


```


* but apparently used, query time cut in half

```
agens=#  MATCH (n:Expression  ) where n.values['ENSG00000000003'] > 50  return count(*) ;
 count
-------
  1357
(1 row)

Time: 30459.542 ms
```



## Traversals

```

MATCH (i:Individual)-[:InProject]->(p:Project) return count(i) ;
 count
-------
   752
(1 row)

Time: 1.990 ms


MATCH (s:BioSample)-[:BiosampleFor]->(i:Individual)-[:InProject]->(p:Project) return count(s) ;
 count
-------
 15598
(1 row)

Time: 1038.061 ms


MATCH (a:Aliquot)-[:AliquotFor]->(s:BioSample)-[:BiosampleFor]->(i:Individual)-[:InProject]->(p:Project) return count(a) ;
 count
-------
 15598
(1 row)

Time: 217.253 ms



explain MATCH (e:Expression)-[:ExpressionOf]->(a:Aliquot)-[:AliquotFor]->(s:BioSample)-[:BiosampleFor]->(i:Individual)-[:InProject]->(p:Project)  return count(e) ;
                                                                                                                                                                           QUERY PLAN

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------
 Aggregate  (cost=7193.74..7193.75 rows=1 width=8)
   ->  Hash Join  (cost=2497.34..7164.53 rows=11684 width=40)
         Hash Cond: ("<0000000013>"."end" = a.id)
         ->  Hash Join  (cost=1897.39..6403.92 rows=11684 width=56)
               Hash Cond: ("<0000000013>".start = e.id)
               ->  Hash Join  (cost=1548.41..5894.29 rows=11684 width=24)
                     Hash Cond: ("<0000000015>"."end" = i.id)
                     Join Filter: ((("<0000000013>".tableoid <> "<0000000016>".tableoid) OR ("<0000000013>".ctid <> "<0000000016>".ctid)) AND (("<0000000014>".tableoid <> "<0000000016>".tableoid) OR ("<0000000014>".ctid <> "<0000000016>".ctid)) AND (("<0000000015>".tableoid <> "<00000
00016>".tableoid) OR ("<0000000015>".ctid <> "<0000000016>".ctid)))
                     ->  Hash Join  (cost=1465.82..5475.72 rows=11686 width=62)
                           Hash Cond: (s.id = "<0000000015>".start)
                           Join Filter: ((("<0000000013>".tableoid <> "<0000000015>".tableoid) OR ("<0000000013>".ctid <> "<0000000015>".ctid)) AND (("<0000000014>".tableoid <> "<0000000015>".tableoid) OR ("<0000000014>".ctid <> "<0000000015>".ctid)))
                           ->  Hash Join  (cost=999.86..4732.20 rows=11687 width=60)
                                 Hash Cond: (s.id = "<0000000014>"."end")
                                 ->  Seq Scan on biosample s  (cost=0.00..3556.98 rows=15598 width=8)
                                 ->  Hash  (cost=853.77..853.77 rows=11687 width=52)
                                       ->  Hash Join  (cost=348.98..853.77 rows=11687 width=52)
                                             Hash Cond: ("<0000000014>".start = "<0000000013>"."end")
                                             Join Filter: (("<0000000013>".tableoid <> "<0000000014>".tableoid) OR ("<0000000013>".ctid <> "<0000000014>".ctid))
                                             ->  Seq Scan on aliquotfor "<0000000014>"  (cost=0.00..270.98 rows=15598 width=26)
                                             ->  Hash  (cost=202.88..202.88 rows=11688 width=26)
                                                   ->  Seq Scan on expressionof "<0000000013>"  (cost=0.00..202.88 rows=11688 width=26)
                           ->  Hash  (cost=270.98..270.98 rows=15598 width=26)
                                 ->  Seq Scan on biosamplefor "<0000000015>"  (cost=0.00..270.98 rows=15598 width=26)
                     ->  Hash  (cost=73.19..73.19 rows=752 width=26)
                           ->  Nested Loop  (cost=22.92..73.19 rows=752 width=26)
                                 Join Filter: ("<0000000016>"."end" = p.id)
                                 ->  Seq Scan on project p  (cost=0.00..1.01 rows=1 width=8)
                                 ->  Hash Join  (cost=22.92..62.78 rows=752 width=34)
                                       Hash Cond: (i.id = "<0000000016>".start)
                                       ->  Seq Scan on individual i  (cost=0.00..29.52 rows=752 width=8)
                                       ->  Hash  (cost=13.52..13.52 rows=752 width=26)
                                             ->  Seq Scan on inproject "<0000000016>"  (cost=0.00..13.52 rows=752 width=26)
               ->  Hash  (cost=202.88..202.88 rows=11688 width=40)
                     ->  Seq Scan on expression e  (cost=0.00..202.88 rows=11688 width=40)
         ->  Hash  (cost=404.98..404.98 rows=15598 width=8)
               ->  Seq Scan on aliquot a  (cost=0.00..404.98 rows=15598 width=8)
(36 rows)



MATCH (e:Expression)-[:ExpressionOf]->(a:Aliquot)-[:AliquotFor]->(s:BioSample)-[:BiosampleFor]->(i:Individual)-[:InProject]->(p:Project)  return count(e) ;
 count
-------
 11688
(1 row)

Time: 32738.211 ms


MATCH (e:Expression)-[:ExpressionOf]->(a:Aliquot)-[:AliquotFor]->(s:BioSample)-[:BiosampleFor]->(i:Individual)-[:InProject]->(p:Project)  return count(e) ;

```


##

```
MATCH
  (e:Expression)-[]->(a:Aliquot)-[]->(s:BioSample)-[]->(i:Individual)-[]->(p:Project)
return
  count(distinct e.gid) as expression,
  count(distinct a.gid) as aliquot,
  count(distinct s.gid) as biosample,
  count(distinct i.gid) as individual,
  count(distinct p.gid) as project
  ;


  MATCH
    (a:Aliquot)-[]->(s:BioSample)-[]->(i:Individual)-[]->(p:Project)
  return
    count(distinct a.gid) as aliquot,
    count(distinct s.gid) as biosample,
    count(distinct i.gid) as individual,
    count(distinct p.gid) as project
    ;


    MATCH
      (i:Individual)-[]->(p:Project)
    return
      count(distinct i.gid) as individual,
      count(distinct p.gid) as project
      ;
```

## unwind not supported https://github.com/bitnine-oss/agensgraph/issues/368
```
MATCH
  (p:Project {gid: 'Project:gtex'}),
  (e:Expression {gid: 'Expression:gtex:GTEX-11GSP-1326-SM-5A5KY'}),
  path = shortestPath((p)-[*..15]-(e))
UNWIND nodes(path) as n
  RETURN { id : id(n), labels : labels(n)} as node, relationships(path) as edges

ERROR:  syntax error at or near "UNWIND"
LINE 5: UNWIND nodes(path) as n

--- workaround

MATCH
  (p:Project {gid: 'Project:gtex'}),
  (e:Expression {gid: 'Expression:gtex:GTEX-11GSP-1326-SM-5A5KY'}),
  path = shortestPath((p)-[*..15]-(e))
  RETURN
  jsonb_pretty(to_jsonb(nodes(path))) as path

  path
--------------------------------------------------------------------------------------------------
[                                                                                               +
{                                                                                           +
"id": "6.1",                                                                            +
"properties": {                                                                         +
"gid": "Project:gtex",                                                              +
"project_id": "gtex",                                                               +
"gdc_attributes": {                                                                 +
}                                                                                   +
}                                                                                       +
},                                                                                          +
{                                                                                           +
"id": "5.11",                                                                           +
"properties": {                                                                         +
"gid": "Individual:GTEX-11GSP",                                                     +
"individual_id": "GTEX-11GSP",                                                      +
"gdc_attributes": {                                                                 +
},                                                                                  +
"gtex_attributes": {                                                                +
"AGE": "60-69",                                                                 +
"SEX": "2",                                                                     +
"SUBJID": "GTEX-11GSP",                                                         +
"DTHHRDY": "2"                                                                  +
}                                                                                   +
}                                                                                       +
},                                                                                          +
{                                                                                           +
"id": "4.574",                                                                          +
"properties": {                                                                         +
"gid": "Biosample:GTEX-11GSP-1326-SM-5A5KY",                                        +
"biosample_id": "GTEX-11GSP-1326-SM-5A5KY",                                         +
"gtex_attributes": {                                                                +
"SMTS": "Heart",                                                                +
"SMGTC": "",                                                                    +


```
