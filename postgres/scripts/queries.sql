

select  distinct label from edge ;


select a.from from edge as a join edge b on (a.from = b.to and a.label = 'AliquotFor') ;

select a.from as "aliquot", a.to as "biosample" from edge as a
where a.label = 'AliquotFor' ;


select a.from as "aliquot", a.to as "biosample", i.to as "individual" , p.to as "project"
  from edge as a
    join edge i on (i.from = a.to and i.label = 'BiosampleFor')
      join edge p on (p.from = i.to and p.label = 'InProject')
  where a.label = 'AliquotFor' ;


select a.from as "aliquot", a.to as "biosample", i.to as "individual" , p.to as "project"
  from edge as a
    join edge i on (a.label = 'AliquotFor' and i.from = a.to and i.label = 'BiosampleFor')
      join edge p on (p.label = 'InProject' and p.from = i.to)
;


select count(distinct a.from) as "aliquot", count(distinct a.to) as "biosample", count(distinct i.to) as "individual" , count(distinct p.to) as "project"
  from edge as a
    join edge i on (a.label = 'AliquotFor' and i.from = a.to and i.label = 'BiosampleFor')
      join edge p on (p.label = 'InProject' and p.from = i.to)
;


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




select
  v.gid as "vertex_gid",
  e.from as "expression",
  a.from as "aliquot",
  a.to as "biosample",
  i.to as "individual",
  p.to as "project"
  from edge as a
    join edge i on (a.label = 'AliquotFor' and i.from = a.to and i.label = 'BiosampleFor')
      join edge p on (p.label = 'InProject' and p.from = i.to)
      join edge e on (e.label = 'ExpressionOf' and a.from = e.to)
      join vertex v on (v.label = 'Expression' and v.gid = e.from)
;
vertex_gid                   |                   expression                   |                aliquot                 |                biosample                 |      individual       |   p
roject
------------------------------------------------+------------------------------------------------+----------------------------------------+------------------------------------------+-----------------------+----
----------
Expression:gtex:GTEX-1117F-0226-SM-5GZZ7       | Expression:gtex:GTEX-1117F-0226-SM-5GZZ7       | Aliquot:GTEX-1117F-0226-SM-5GZZ7       | Biosample:GTEX-1117F-0226-SM-5GZZ7       | Individual:GTEX-1117F | Pro
ject:gtex
Expression:gtex:GTEX-1117F-0426-SM-5EGHI       | Expression:gtex:GTEX-1117F-0426-SM-5EGHI       | Aliquot:GTEX-1117F-0426-SM-5EGHI       | Biosample:GTEX-1117F-0426-SM-5EGHI       | Individual:GTEX-1117F | Pro
ject:gtex
Expression:gtex:GTEX-1117F-0526-SM-5EGHJ       | Expression:gtex:GTEX-1117F-0526-SM-5EGHJ       | Aliquot:GTEX-1117F-0526-SM-5EGHJ       | Biosample:GTEX-1117F-0526-SM-5EGHJ       | Individual:GTEX-1117F | Pro
ject:gtex....
Time: 325.977 ms



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
;

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


select
  v.gid as "vertex_gid",
  v.data->'metric' as "metric",
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
;
Time: 37637.301 ms


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
      join vertex v on (v.gid = e.from)
where
  cast (v.data->'values'->>'ENSG00000227232' as float ) > 15 ;
  and v.label = 'Expression'
;


select
  v.gid as "vertex_gid",
  v.data->'values'->'ENSG00000227232' as "ENSG00000227232",
  e.to
  from vertex as v
    join edge e on (e.label = 'ExpressionOf' and v.gid = e.from)
where
  v.label = 'Expression'
  and
  cast (v.data->'values'->>'ENSG00000227232' as float ) > 15
;

select
  v.gid as "vertex_gid",
  v.data->'values'->'ENSG00000227232' as "ENSG00000227232"
  from vertex as v
where
  v.label = 'Expression'
  and
  cast (v.data->'values'->>'ENSG00000227232' as float ) > 15
;
Time: 15608.100 ms

select
  v.gid as "vertex_gid",
  v.data->'values'->>'ENSG00000227232' as "ENSG00000227232",
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
;
Time: 44902.668 ms


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
;
46150.486 ms



select
  gid,
  data->'values'->>'ENSG00000227232' as "ENSG00000227232"
from vertex as v
  where v.label = 'Expression'
  and (v.data->'values'->>'ENSG00000227232')::float > '15'  limit 1;

select
  avg((v.data->'values'->>'ENSG00000227232')::float) as "average_ENSG00000227232"
from vertex as v
  where v.label = 'Expression' ;


CREATE INDEX vertex_ENSG00000227232 ON vertex (cast (data->'values'->>'ENSG00000227232' as float) );


# select gid, v.data->'values'->'ENSG00000227232' as "ENSG00000227232"
from vertex as v
 where v.label = 'Expression'
  and v.data->'values'->'ENSG00000227232' > '15'::jsonb ;
  gid                       |  ENSG00000227232
------------------------------------------------+--------------------
Expression:gtex:GTEX-15DCZ-1726-SM-69LPX       | 15.630000114440918
Expression:gtex:GTEX-15G19-1126-SM-7KUEL       | 16.770000457763672
Expression:gtex:GTEX-15SHW-1726-SM-6LPJ2       | 16.309999465942383
Expression:gtex:GTEX-16NPV-0626-SM-7IGOR       | 22.56999969482422
Expression:gtex:GTEX-17F96-0726-SM-793CB       | 15.600000381469727
...
Time: 15439.542 ms


# select avg((v.data->'values'->>'ENSG00000227232')::float) as "ENSG00000227232"
from vertex as v
where v.label = 'Expression'
and v.data->'values'->'ENSG00000227232' > '15'::jsonb ;

   ENSG00000227232
  ------------------
   21.1918261183091
  (1 row)

  Time: 14799.759 ms




# CREATE INDEX vertex_data on vertex USING gin (data);
CREATE INDEX
Time: 5797025.542 ms





CREATE INDEX vertex_gtex_values ON
  vertex (
    cast (data->'values'->>'value' as float),
    cast (data->'values'->>'name' as varchar)
  );


select
  gid,
  data->'values'->1->>'name'
from vertex as v
  where v.label = 'Expression'
limit 1 ;


select
  gid,
  jsonb_array_elements(data->'values')->>'name'
from vertex as v
  where v.label = 'Expression'
    and v.data->'values'->'value' > '15'::jsonb
limit 5 ;


explain select
  gid,
  jsonb_array_elements(data->'values')->>'name',
  jsonb_array_elements(data->'values')->>'value'
from vertex as v
  where v.label = 'Expression'
    and cast (v.data->'values'->>'value' as float ) > 15 ;

                                                QUERY PLAN
-----------------------------------------------------------------------------------------------------------
 Index Scan using vertex_gtex_values on vertex v  (cost=0.29..8.82 rows=100 width=99)
   Index Cond: ((((data -> 'values'::text) ->> 'value'::text))::double precision > '15'::double precision)
   Filter: ((label)::text = 'Expression'::text)
(3 rows)

select
  gid,
  jsonb_array_elements(data->'values')->>'name',
  jsonb_array_elements(data->'values')->>'value'
from vertex as v
  where v.label = 'Expression'
    and cast (v.data->'values'->>'value' as float ) > 9
limit 10
    ;


http://erthalion.info/2017/12/21/advanced-json-benchmarks/
http://niquola.github.io/blog/postgresql/2014/08/20/pg-jsquery-2014.html

select
  gid,
  jsonb_array_elements(data->'values')->>'name',
  jsonb_array_elements(data->'values')->>'value',
  cast (v.data->'values'->>'value' as float )
from vertex as v
  where v.label = 'Expression'
limit 10
    ;


show sizes

```

    CREATE OR REPLACE FUNCTION bmeg() RETURNS TABLE(relation varchar, kind varchar, size varchar)
        AS $$
        SELECT nspname || '.' || relname AS "relation",
            relkind as "kind",
            pg_size_pretty(pg_relation_size(C.oid)) AS "size"
          FROM pg_class C
          LEFT JOIN pg_namespace N ON (N.oid = C.relnamespace)
          WHERE nspname NOT IN ('pg_catalog', 'information_schema',  'pg_toast')
          ORDER BY pg_relation_size(C.oid) DESC
          LIMIT 20
        $$
        LANGUAGE SQL;


        select * from bmeg();


        bmeg_test=# select * from bmeg();
           relation    |  size
        ---------------+--------
         public.vertex | 832 MB
         public.edge   | 753 MB
        (2 rows)
```



show duplicates

```

select * from vertex out
where (select count(*) from vertex as inn
where inn.gid = out.gid) > 1
limit 1 ;
```


montor index build
```
SELECT
  t.tablename,
  indexname,
  c.reltuples::int AS num_rows,
  pg_size_pretty(pg_relation_size(quote_ident(t.tablename)::text)) AS table_size,
  pg_size_pretty(pg_relation_size(quote_ident(indexrelname)::text)) AS index_size,
  CASE WHEN indisunique THEN 'Y'
    ELSE 'N'
  END AS UNIQUE,
  idx_scan AS number_of_scans,
  idx_tup_read AS tuples_read,
  idx_tup_fetch AS tuples_fetched
FROM pg_tables t
  LEFT OUTER JOIN pg_class c ON t.tablename=c.relname
  LEFT OUTER JOIN
    ( SELECT c.relname AS ctablename, ipg.relname AS indexname, x.indnatts AS number_of_columns, idx_scan, idx_tup_read, idx_tup_fetch, indexrelname, indisunique FROM pg_index x
      JOIN pg_class c ON c.oid = x.indrelid
      JOIN pg_class ipg ON ipg.oid = x.indexrelid
      JOIN pg_stat_all_indexes psai ON x.indexrelid = psai.indexrelid )
    AS foo
  ON t.tablename = foo.ctablename
WHERE t.schemaname='public'
ORDER BY 1,2;
```
