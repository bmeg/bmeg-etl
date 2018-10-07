

# BMEG SQL Data Model and Implementation

## Scope

* Evaluate postgres support for BMEG's graph oriented `vertex` and `edge` types.
* Evaluate postgres support for BMEG's document model, specifically json queries

## Out of Scope

* "Bake-off" with mongo
* Translation and strategies for "graph to sql" i.e grip queries query translation (all queries were hand crafted)

## Results

* Load into postgres trivial
* Graph indexing and traversal queries performant
* Querying vertex document complex and performance inconclusive

## Details

### Data Model

The postgres model mirrors the mongo implementation, two simple tables, sharded by [label, gid] each with a data column that contains a json document with all entity properties.

![image](https://user-images.githubusercontent.com/47808/45767329-bbaba900-bbee-11e8-85ae-a6f764327fb5.png)


```
CREATE TABLE IF NOT EXISTS  vertex (
 gid varchar not null,
 label varchar not null,
 data jsonb
);

CREATE TABLE IF NOT EXISTS  edge (
 gid varchar not null,
 label varchar not null,
 "from" varchar not null,
 "to" varchar not null,
 data jsonb
);
```

### Load & Indexing

A simple python [script](https://github.com/bmeg/bmeg-etl/blob/postgres/postgres/scripts/loader.py) leverages [sqlalchemy](https://www.sqlalchemy.org/) and [dataset](https://dataset.readthedocs.io/en/latest/quickstart.html)  to load the tables.


The graph traversal indexes are likewise simple and straightforward.

```
CREATE INDEX vertex_gid ON vertex (gid);
CREATE INDEX vertex_label ON vertex (label);
CREATE INDEX edge_label_from_to ON edge (label, "from", "to");
CREATE INDEX edge_label_to_from ON edge (label, "to", "from");
```

### Test Data

The GTEx data set was loaded.  It was selected because:
* shares a common `Project->BioSample->Aliquot->...` model with many other data sources
* has a heavyweight `Expresssion` vertex with a 52K column gene expression

### Graph traversal queries

```
select
  count(distinct e.to) as "expression",
  count(distinct a.from) as "aliquot",
  count(distinct a.to) as "biosample",
  count(distinct i.to) as "individual",
  count(distinct p.to) as "project"
  from edge as a
    join edge i on (a.label = 'AliquotFor' and i.from = a.to and i.label = 'BiosampleFor')
      join edge p on (p.label = 'InProject' and p.from = i.to)
      join edge e on (e.label = 'ExpressionOf' and a.from = e.to)
;

expression | aliquot | biosample | individual | project
------------+---------+-----------+------------+---------
     11688 |   11688 |     11688 |        714 |       1
(1 row)

Time: 169.941 ms
```

### Vertex data jsonb queries

Postgres operators are adequate for simple queries.
Where:
* `->` return the jsonb type
* `->>` return as sql text


```
select
  v.gid as "vertex_gid",
  v.data->'gtex_attributes'->>'AGE' as "age",
  e.from as "expression",
  a.from as "aliquot",
  a.to as "biosample",
  i.to as "individual",
  p.to as "project"
  from edge as a
    join edge i on (a.label = 'AliquotFor' and i.from = a.to and i.label = 'BiosampleFor')
      join edge p on (p.label = 'InProject' and p.from = i.to)
      join edge e on (e.label = 'ExpressionOf' and a.from = e.to)
      join vertex v on (v.label = 'Individual' and v.gid = i.to)
  where v.data->'gtex_attributes'->>'AGE' = '60-69'
  and a.from = 'Aliquot:GTEX-131XF-1426-SM-5BC68'
;
```

It get more complicated for numeric values.
i.e.  `::jsonb numeric` is not the same type as sql `float` and cannot be used without casting


```
select
  v.gid as "vertex_gid",
  v.data->'values'->'ENSG00000227232' as "ENSG00000227232",
  e.to as "expression",
  a.from as "aliquot",
  a.to as "biosample",
  i.to as "individual",
  p.to as "project"
  from edge as a
    join edge i on (a.label = 'AliquotFor' and i.from = a.to and i.label = 'BiosampleFor')
      join edge p on (p.label = 'InProject' and p.from = i.to)
      join edge e on (e.label = 'ExpressionOf' and a.from = e.to)
      join vertex v on (v.label = 'Expression' and v.gid = e.from)
where
  v.data->'values'->'ENSG00000227232' > '15'::jsonb ;
;
Time: 64801.650 ms

select
  v.gid as "vertex_gid",
  v.data->'values'->'ENSG00000227232' as "ENSG00000227232",
  e.to as "expression",
  a.from as "aliquot",
  a.to as "biosample",
  i.to as "individual",
  p.to as "project"
  from edge as a
    join edge i on (a.label = 'AliquotFor' and i.from = a.to and i.label = 'BiosampleFor')
      join edge p on (p.label = 'InProject' and p.from = i.to)
      join edge e on (e.label = 'ExpressionOf' and a.from = e.to)
      join vertex v on (v.label = 'Expression' and v.gid = e.from)
where
  cast (v.data->'values'->>'ENSG00000227232' as float ) > 15 ;
;
Time: 61650.050 ms


select avg((v.data->'values'->>'ENSG00000227232')::float) as "ENSG00000227232"
from vertex as v
where v.label = 'Expression'
and v.data->'values'->'ENSG00000227232' > '15'::jsonb ;

   ENSG00000227232
  ------------------
   21.1918261183091
  (1 row)

  Time: 14799.759 ms

```
Note:  this will get easier with postgres 11 (currently in beta)

![image](https://user-images.githubusercontent.com/47808/45837029-33e39e80-bcc3-11e8-89c1-d08090b49548.png)




### jsonb Indexing

GIN indexes can be used to efficiently search for keys or key/value pairs.

```
CREATE INDEX vertex_data on vertex USING gin (data);
```

However, they do not seem to be useful at all for numeric comparisons.  Checks only equality for
scalars.

This requires additional ad-hoc functional indexes.  e.g.

```
# return data.values.ENSG00000227232 as text and cast it to float
CREATE INDEX vertex_ENSG00000227232 ON vertex (cast (data->'values'->>'ENSG00000227232' as float) );
```

### JSQuery

Stock postgres json operators do not enable access to json arrays.  

Fortunately, a postgres plugin [jsquery](http://www.sai.msu.su/~megera/postgres/talks/pgconfeu-2014-jsquery.pdf) enables a complete set of jsonb operators.  A [pre-built binary](https://packages.ubuntu.com/bionic/postgresql-10-jsquery) for ubuntu exists.

![image](https://user-images.githubusercontent.com/47808/45771018-d97e0b80-bbf8-11e8-9e62-8ab85fad3292.png)

![image](https://user-images.githubusercontent.com/47808/45771098-0d593100-bbf9-11e8-9c41-48845e88e89b.png)




Postgres json query performance on unindexed field fairly poor.
```
test=# explain
select gid, v.data->'values'->'ENSG00000227232' as "ENSG00000227232"
from vertex as v
where v.label = 'Expression'
and v.data->'values'->'ENSG00000227232' > '15'::jsonb ;
                                   QUERY PLAN
---------------------------------------------------------------------------------
 Bitmap Heap Scan on vertex v  (cost=257.20..4234.24 rows=3908 width=67)
   Recheck Cond: ((label)::text = 'Expression'::text)
   Filter: (((data -> 'values'::text) -> 'ENSG00000227232'::text) > '15'::jsonb)
   ->  Bitmap Index Scan on vertex_label  (cost=0.00..256.23 rows=11725 width=0)
         Index Cond: ((label)::text = 'Expression'::text)
(5 rows)

...
Time: 90621.367 ms (01:30.621)

```

JSQuery performance on unindexed field much better.
```
test=# explain
select gid, v.data->'values'->'ENSG00000227232' as "ENSG00000227232"
from vertex as v
where v.label = 'Expression'
and
v.data @@ 'values.ENSG00000227232 > 15' ;

QUERY PLAN
---------------------------------------------------------------------------------
Bitmap Heap Scan on vertex v  (cost=256.23..4155.17 rows=12 width=67)
Recheck Cond: ((label)::text = 'Expression'::text)
Filter: (data @@ '"values"."ENSG00000227232" > 15'::jsquery)
->  Bitmap Index Scan on vertex_label  (cost=0.00..256.23 rows=11725 width=0)
Index Cond: ((label)::text = 'Expression'::text)
(5 rows)

...
Time: 44807.943 ms (00:44.808)

```

Create a index on single column using ::jsonb


`CREATE INDEX vertex_ENSG00000227232 ON vertex ((data->'values'->'ENSG00000227232'));`

Index not used by jsquery

```
explain
select gid,  v.data->'values'->'ENSG00000227232' as "ENSG00000227232"
from vertex as v
where
v.data @@ 'values.ENSG00000227232 > 15'
and
v.label = 'Expression' ;

QUERY PLAN
---------------------------------------------------------------------------------
Bitmap Heap Scan on vertex v  (cost=256.43..4155.72 rows=12 width=36)
Recheck Cond: ((label)::text = 'Expression'::text)
Filter: (data @@ '"values"."ENSG00000227232" > 15'::jsquery)
->  Bitmap Index Scan on vertex_label  (cost=0.00..256.42 rows=11751 width=0)
Index Cond: ((label)::text = 'Expression'::text)
(5 rows
...
38416.829 ms (00:38.417)

```


Used by postgres json

```
explain
select gid, v.data->'values'->'ENSG00000227232' as "ENSG00000227232"
from vertex as v
where v.label = 'Expression'
and v.data->'values'->'ENSG00000227232' > '15'::jsonb ;
QUERY PLAN
----------------------------------------------------------------------------------------------------------------------------------
Bitmap Heap Scan on vertex v  (cost=341.52..2617.64 rows=1007 width=67)
Recheck Cond: ((((data -> 'values'::text) -> 'ENSG00000227232'::text) > '15'::jsonb) AND ((label)::text = 'Expression'::text))
->  BitmapAnd  (cost=341.52..341.52 rows=1007 width=0)
->  Bitmap Index Scan on vertex_ensg00000227232  (cost=0.00..84.34 rows=3740 width=0)
Index Cond: (((data -> 'values'::text) -> 'ENSG00000227232'::text) > '15'::jsonb)
->  Bitmap Index Scan on vertex_label  (cost=0.00..256.42 rows=11751 width=0)
Index Cond: ((label)::text = 'Expression'::text)''

...

Time: 9906.712 ms (00:09.907)

```


Best performing query


Create native float index
`CREATE INDEX vertex_ENSG00000227232 ON vertex (cast (data->'values'->>'ENSG00000227232' as float) );`


```
explain
select gid,  v.data->'values'->'ENSG00000227232' as "ENSG00000227232"
from vertex as v
where
cast (data->'values'->>'ENSG00000227232' as float) > 15
and
v.label = 'Expression' ;

QUERY PLAN
------------------------------------------------------------------------------------------------------------------------------------------------------------------
Bitmap Heap Scan on vertex v  (cost=323.28..2533.61 rows=959 width=67)
Recheck Cond: (((((data -> 'values'::text) ->> 'ENSG00000227232'::text))::double precision > '15'::double precision) AND ((label)::text = 'Expression'::text))
->  BitmapAnd  (cost=323.28..323.28 rows=959 width=0)
->  Bitmap Index Scan on vertex_ensg00000227232  (cost=0.00..71.38 rows=3612 width=0)
Index Cond: ((((data -> 'values'::text) ->> 'ENSG00000227232'::text))::double precision > '15'::double precision)
->  Bitmap Index Scan on vertex_label  (cost=0.00..251.17 rows=11584 width=0)
Index Cond: ((label)::text = 'Expression'::text)
...
Time: 10508.383 ms (00:10.508)


```
