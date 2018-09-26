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

* not shown in explain ?

```
 agens=# explain MATCH (n:Expression  ) where n.values['ENSG00000000003'] > 50  return count(*) ;
                                       QUERY PLAN
 ---------------------------------------------------------------------------------------
  Aggregate  (cost=241.84..241.85 rows=1 width=8)
    ->  Seq Scan on expression n  (cost=0.00..232.10 rows=3896 width=0)
          Filter: (properties.'values'::text['"ENSG00000000003"'::jsonb] > '50'::jsonb)
 (3 rows)

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
