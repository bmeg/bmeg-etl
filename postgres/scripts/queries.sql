

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


turn off paging
```
\pset pager off
```


show sizes

```

CREATE OR REPLACE FUNCTION sizes() RETURNS TABLE(relation varchar, kind varchar, size varchar)
  AS $$
  SELECT nspname || '.' || relname AS "relation",
      relkind || '' as "kind",
      pg_size_pretty(pg_relation_size(C.oid)) AS "size"
    FROM pg_class C
    LEFT JOIN pg_namespace N ON (N.oid = C.relnamespace)
    WHERE nspname NOT IN ('pg_catalog', 'information_schema',  'pg_toast')
    ORDER BY pg_relation_size(C.oid) DESC
    LIMIT 20
  $$
  LANGUAGE SQL;


select * from sizes();


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


Allele:ab5f4bf29886c718f4b5ed285c7bbf838a0e91d4


select
  count(distinct a.from) as "aliquot",
  count(distinct a.to) as "biosample",
  count(distinct i.to) as "individual",
  count(distinct p.to) as "project"
  from edge as a
    join edge i on (a.label = 'AliquotFor' and i.from = a.to and i.label = 'BiosampleFor')
      join edge p on (p.label = 'InProject' and p.from = i.to)
;



select
  count(distinct al.from) as "allele",
  count(distinct a.from) as "aliquot",
  count(distinct a.to) as "biosample",
  count(distinct i.to) as "individual",
  count(distinct p.to) as "project"
  from edge as a
    join edge i on (a.label = 'AliquotFor' and i.from = a.to and i.label = 'BiosampleFor')
      join edge p on (p.label = 'InProject' and p.from = i.to)
      join edge al on (al.label = 'AlleleCall' and a.from = al.to)
;



select
  count(distinct al.from) as "allele",
  count(distinct c.to) as "callset",
  count(distinct a.from) as "aliquot",
  count(distinct a.to) as "biosample",
  count(distinct i.to) as "individual",
  count(distinct p.to) as "project"
  from edge as a
    join edge i on (a.label = 'AliquotFor' and i.from = a.to and i.label = 'BiosampleFor')
      join edge p on (p.label = 'InProject' and p.from = i.to)
        join edge c on (c.label = 'CallsetFor' and a.from = c.to and a.label = 'AliquotFor')
          join edge al on (al.label = 'AlleleCall' and c.from = al.to and c.label = 'CallsetFor')
;




explain select
  count(distinct al.from) as "allele",
  count(distinct c.to) as "callset",
  count(distinct a.from) as "aliquot",
  count(distinct a.to) as "biosample",
  count(distinct i.to) as "individual",
  count(distinct p.to) as "project"
  from edge as a
    join edge i on (a.label = 'AliquotFor' and i.from = a.to and i.label = 'BiosampleFor')
      join edge p on (p.label = 'InProject' and p.from = i.to)
        join edge c on (c.label = 'CallsetFor' and a.from = c.to)
          join edge al on (al.label = 'AlleleCall' and c.from = al.to)
;



explain select
  count(distinct al.from) as "allele",
  count(distinct c.to) as "callset",
  count(distinct a.from) as "aliquot",
  count(distinct a.to) as "biosample",
  count(distinct i.to) as "individual",
  count(distinct p.to) as "project"
  from edge as a
    join edge i on (a.label = 'AliquotFor' and i.from = a.to and i.label = 'BiosampleFor')
      join edge p on (p.label = 'InProject' and p.from = i.to)
        join edge c on (c.label = 'CallsetFor' and a.from = c.to)
          join edge al on (al.label = 'AlleleCall' and c.from = al.to)
  WHERE
    al.from = 'Allele:ab5f4bf29886c718f4b5ed285c7bbf838a0e91d4'
;

select
  count(distinct al.from) as "allele",
  count(distinct c.to) as "callset",
  count(distinct a.from) as "aliquot",
  count(distinct a.to) as "biosample",
  count(distinct i.to) as "individual",
  p.to as "project"
  from edge as a
    join edge i on (a.label = 'AliquotFor' and i.from = a.to and i.label = 'BiosampleFor')
      join edge p on (p.label = 'InProject' and p.from = i.to)
        join edge c on (c.label = 'CallsetFor' and a.from = c.to)
          join edge al on (al.label = 'AlleleCall' and c.from = al.to)
  WHERE
    al.from = 'Allele:ab5f4bf29886c718f4b5ed285c7bbf838a0e91d4'
  GROUP BY
    project
;


select
  count(distinct al.from) as "allele",
  count(distinct c.to) as "callset",
  count(distinct a.from) as "aliquot",
  count(distinct a.to) as "biosample",
  count(distinct i.to) as "individual",
  p.to as "project"
  from edge as a
    join edge i on (a.label = 'AliquotFor' and i.from = a.to and i.label = 'BiosampleFor')
      join edge p on (p.label = 'InProject' and p.from = i.to)
        join edge c on (c.label = 'CallsetFor' and a.from = c.to)
          join edge al on (al.label = 'AlleleCall' and c.from = al.to)
  GROUP BY
    project
;



Show the most 'popular' alleles

```

select
  al.from  as "allele",
  count(distinct c.to) as "callset",
  count(distinct a.from) as "aliquot",
  count(distinct a.to) as "biosample",
  count(distinct i.to) as "individual",
  count(distinct p.to) as "project"
  from edge as a
    join edge i on (a.label = 'AliquotFor' and i.from = a.to and i.label = 'BiosampleFor')
      join edge p on (p.label = 'InProject' and p.from = i.to)
        join edge c on (c.label = 'CallsetFor' and a.from = c.to)
          join edge al on (al.label = 'AlleleCall' and c.from = al.to)
  GROUP BY
    allele
  HAVING
   count(distinct p.to)  > 10
  ORDER  BY project DESC;

                     allele                      | callset | aliquot | biosample | individual | project
-------------------------------------------------+---------+---------+-----------+------------+---------
 Allele:ebc8cd9acf2198d45dbdfd764a69368dca489a0b |     296 |     296 |       296 |        164 |      37
 Allele:8e9f5f0e7766550b409878711e929445645ba8c8 |     354 |     354 |       354 |        194 |      36
 Allele:d54fe7c2a96ce7ac0f0f75096926c1436d61e5ad |     314 |     314 |       314 |        184 |      35
 Allele:7d9fbcd1c09166fab2ea01d1bdd571255ddd830c |     210 |     210 |       210 |        116 |      34
 Allele:2f29241c6df6db2090639f153379e3da7e1cf8e9 |     569 |     569 |       569 |        292 |      33
 Allele:74db52a85f0d48bf3b1e6e6ebb10c68d9e929ca0 |     257 |     257 |       257 |        148 |      32
 Allele:2d5c0debfaf3ba7a5e5e0614203512ff0a3a52d5 |     176 |     176 |       176 |         94 |      32
 Allele:dd9394f9294c03b23da21f612bb88f22f40b6f87 |     469 |     469 |       469 |        263 |      30
 Allele:ce45202f2eb047f1b22c13fb17c5f84cbb06e195 |     531 |     531 |       531 |        273 |      29
 Allele:18f46c782a347c0ed6f7c454fcd7ad9c5e209965 |     160 |     160 |       160 |         89 |      29
 Allele:e0cd8bc1a8be06e6a8d44c5fda6748de27b7b666 |     107 |     107 |       107 |         62 |      26
 Allele:9a7be53d7094bef75f035c264a1f4e25896b440e |     173 |     173 |       173 |        104 |      25
 Allele:4fc7d3a2bde4e19cd8c6260adc7ec8a252898eef |     395 |     395 |       395 |        213 |      24
 Allele:26b1a6cb66bc58350ca21ef592b26d075a0a7a0f |     243 |     243 |       243 |        131 |      23
 Allele:b52ced1e915b178dd364afac3458a17f6e65630c |     358 |     358 |       358 |        183 |      23
 Allele:7fd5a7a99a15dd9695f369b656d15f008f9dad58 |      79 |      79 |        79 |         71 |      23
 Allele:1fd3ee6ba5430e1702483055b912a6bda6db87c8 |     187 |     187 |       187 |         96 |      23
 Allele:6794c0e81a961346669f667e696e566c7767059f |    1248 |    1248 |      1248 |        660 |      23
 Allele:5d7a7f15245108d1c93e06df4442a8b6eac209ec |     254 |     254 |       254 |        158 |      22
 Allele:8bdb32c48766f14bd6055150be4b6b1c2faf401a |      75 |      75 |        75 |         41 |      22
 Allele:17a1c142dbae32abb5c5b11893aaeb31cd602c3d |     244 |     244 |       244 |        125 |      21
 Allele:64882d92d160f0c43ad9cd427745fa83849777c6 |      88 |      88 |        88 |         47 |      21
 Allele:494394b3839557c97fd380f96e94b39c23e4072e |      56 |      56 |        56 |         56 |      21
 Allele:c87423aed4064e419cee7eb1d18564474dfadb34 |      53 |      53 |        53 |         51 |      21
 Allele:b4aa2bfe2b2d6fd36d7756995c9024ce0f7881ee |      49 |      49 |        49 |         49 |      21

 ....

 Time: 116398.436 ms (01:56.398)
 ```


 Alleles with g2p entries
 ```
 select
   al.to  as "allele",
   count(distinct al.from) as "g2p_association"
   from edge as al
 WHERE
   al.label = 'HasAlleleFeature'
 GROUP BY
   allele
 HAVING
   count(distinct al.from)  > 10
 ORDER BY
   "g2p_association" DESC
 ;
```


alleles that exist in a project and exist in g2p

```

select
  count(distinct al.from)  as "allele",
  count(distinct c.to) as "callset",
  count(distinct a.from) as "aliquot",
  count(distinct a.to) as "biosample",
  count(distinct i.to) as "individual",
  count(distinct p.to) as "project",
  count(distinct g2p.from) as "g2p_association"
  from edge as a
    join edge i on (a.label = 'AliquotFor' and i.from = a.to and i.label = 'BiosampleFor')
      join edge p on (p.label = 'InProject' and p.from = i.to)
        join edge c on (c.label = 'CallsetFor' and a.from = c.to)
          join edge al on (al.label = 'AlleleCall' and c.from = al.to)
            join edge g2p on (g2p.label = 'HasAlleleFeature' and g2p.to = al.from)
;

allele | callset | aliquot | biosample | individual | project | g2p_association
--------+---------+---------+-----------+------------+---------+-----------------
  1385 |   12280 |   12280 |     12280 |       6523 |      67 |            8314
(1 row)

Time: 44967.470 ms (00:44.967)

```



select
  count(distinct al.from)  as "allele",
  count(distinct c.to) as "callset",
  count(distinct a.from) as "aliquot",
  count(distinct a.to) as "biosample",
  count(distinct i.to) as "individual",
  count(distinct p.to) as "project",
  count(distinct g2p.from) as "g2p_association",
  count(distinct ph.to) as "phenotype"
  from edge as a
    join edge i on (a.label = 'AliquotFor' and i.from = a.to and i.label = 'BiosampleFor')
      join edge p on (p.label = 'InProject' and p.from = i.to)
        join edge c on (c.label = 'CallsetFor' and a.from = c.to)
          join edge al on (al.label = 'AlleleCall' and c.from = al.to)
            join edge g2p on (g2p.label = 'HasAlleleFeature' and g2p.to = al.from)
              join edge ph on (ph.label = 'HasPhenotype' and g2p.from = ph.from)

;

```
select
  count(distinct al.from)  as "allele",
  count(distinct c.to) as "callset",
  count(distinct a.from) as "aliquot",
  count(distinct a.to) as "biosample",
  count(distinct i.to) as "individual",
  count(distinct p.to) as "project",
  count(distinct g2p.from) as "g2p_association",
  count(distinct ph.to) as "phenotype",
  count(distinct e.to) as "environment"
  from edge as a
    join edge i on (a.label = 'AliquotFor' and i.from = a.to and i.label = 'BiosampleFor')
      join edge p on (p.label = 'InProject' and p.from = i.to)
        join edge c on (c.label = 'CallsetFor' and a.from = c.to)
          join edge al on (al.label = 'AlleleCall' and c.from = al.to)
            join edge g2p on (g2p.label = 'HasAlleleFeature' and g2p.to = al.from)
              join edge ph on (ph.label = 'HasPhenotype' and g2p.from = ph.from)
              join edge e on (e.label = 'HasEnvironment' and g2p.from = e.from)
;
allele | callset | aliquot | biosample | individual | project | g2p_association | phenotype | environment
--------+---------+---------+-----------+------------+---------+-----------------+-----------+-------------
   372 |    8651 |    8651 |      8651 |       4608 |      63 |             247 |       314 |         289
(1 row)

Time: 36476.066 ms (00:36.476)
```

what are the counts of all alleles :
* that have individuals in projects who have had treatments
* that are have evidence with a  compound

```
select
  count(distinct al.from)  as "allele",
  count(distinct c.to) as "callset",
  count(distinct a.from) as "aliquot",
  count(distinct a.to) as "biosample",
  count(distinct i.to) as "individual",
  count(distinct p.to) as "project",
  count(distinct g2p.from) as "g2p_association",
  count(distinct ph.to) as "phenotype",
  count(distinct e.to) as "environment",
  count(distinct t.to) as "treatment"
  from edge as a
    join edge i on (a.label = 'AliquotFor' and i.from = a.to and i.label = 'BiosampleFor')
      join edge p on (p.label = 'InProject' and p.from = i.to)
      join edge t on (t.label = 'TreatedWith' and t.from = i.to)
        join edge c on (c.label = 'CallsetFor' and a.from = c.to)
          join edge al on (al.label = 'AlleleCall' and c.from = al.to)
            join edge g2p on (g2p.label = 'HasAlleleFeature' and g2p.to = al.from)
              join edge ph on (ph.label = 'HasPhenotype' and g2p.from = ph.from)
              join edge e on (e.label = 'HasEnvironment' and g2p.from = e.from)
              join vertex g2p_v on (g2p_v.label = 'G2PAssociation' and g2p.from = g2p_v.to)

;

allele | callset | aliquot | biosample | individual | project | g2p_association | phenotype | environment | treatment
--------+---------+---------+-----------+------------+---------+-----------------+-----------+-------------+-----------
   174 |    1878 |    1878 |      1878 |        908 |      19 |             102 |       163 |         115 |       150
(1 row)

Time: 24250.170 ms (00:24.250)

```


set maintenance_work_mem = '10GB' ;
set work_mem = '10GB' ;
set max_parallel_workers_per_gather = 12 ;
set effective_io_concurrency = 12 ;

select
  count(distinct al.from)  as "allele",
  count(distinct c.to) as "callset",
  count(distinct a.from) as "aliquot",
  count(distinct a.to) as "biosample",
  count(distinct i.to) as "individual",
  count(distinct p.to) as "project",
  count(distinct g2p.from) as "g2p_association",
  count(distinct ph.to) as "phenotype",
  count(distinct e.to) as "environment",
  count(distinct t.to) as "treatment",
  g2p_v.data->'source' as "source"
  from edge as a
    join edge i on (a.label = 'AliquotFor' and i.from = a.to and i.label = 'BiosampleFor')
      join edge p on (p.label = 'InProject' and p.from = i.to)
      join edge t on (t.label = 'TreatedWith' and t.from = i.to)
        join edge c on (c.label = 'CallsetFor' and a.from = c.to)
          join edge al on (al.label = 'AlleleCall' and c.from = al.to)
            join edge g2p on (g2p.label = 'HasAlleleFeature' and g2p.to = al.from)
              join edge ph on (ph.label = 'HasPhenotype' and g2p.from = ph.from)
              join edge e on (e.label = 'HasEnvironment' and g2p.from = e.from)
              join vertex g2p_v on (g2p_v.label = 'G2PAssociation' and g2p.from = g2p_v.gid)
  GROUP BY
  g2p_v.data->'source'
;




"""
select
    {SELECT}
  from edge as a
  {WHERE}
  {GROUP_BY}
  {ORDER_BY}
;
""".replace

def Q(q):
    join = """
    edge as a
      join edge i on (a.label = 'AliquotFor' and i.from = a.to and i.label = 'BiosampleFor')
        join edge p on (p.label = 'InProject' and p.from = i.to)
        join edge t on (t.label = 'TreatedWith' and t.from = i.to)
          join edge c on (c.label = 'CallsetFor' and a.from = c.to)
            join edge al on (al.label = 'AlleleCall' and c.from = al.to)
              join edge g2p on (g2p.label = 'HasAlleleFeature' and g2p.to = al.from)
                join edge ph on (ph.label = 'HasPhenotype' and g2p.from = ph.from)
                join edge e on (e.label = 'HasEnvironment' and g2p.from = e.from)
                join vertex g2p_v on (g2p_v.label = 'G2PAssociation' and g2p.from = g2p_v.gid)
    """
    return q.format(join=join)


view of raw edge joins

```
create or replace view alleles
  as
  select
    p.to      as "project_gid",
    i.to      as "individual_gid",
    i.from    as "biosample_gid",
    c.to      as "aliquot_gid",
    c.from    as "callset_gid",
    al.from   as "allele_gid",
    g2p.from  as "g2passociation_gid",
    ph.to     as "phenotype_gid",
    e.to      as "environment_gid",
    t.to      as "treatment_gid"
  from
    edge as a
      join edge i on (a.label = 'AliquotFor' and i.from = a.to and i.label = 'BiosampleFor')
        join edge p on (p.label = 'InProject' and p.from = i.to)
        join edge t on (t.label = 'TreatedWith' and t.from = i.to)
          join edge c on (c.label = 'CallsetFor' and a.from = c.to)
            join edge al on (al.label = 'AlleleCall' and c.from = al.to)
              join edge g2p on (g2p.label = 'HasAlleleFeature' and g2p.to = al.from)
                join edge ph on (ph.label = 'HasPhenotype' and g2p.from = ph.from)
                join edge e on (e.label = 'HasEnvironment' and g2p.from = e.from)
                join vertex g2p_v on (g2p_v.label = 'G2PAssociation' and g2p.from = g2p_v.gid)
;

```

use the edge view to conviently write queries

```
select
  g2p.data->>'source' as "source",
  g2p.data->>'evidence_label' as "evidence_label",
  count(distinct alleles.allele_gid) as "allele_count"
from
  alleles
    join vertex g2p on (g2p.label = 'G2PAssociation' and alleles.g2passociation_gid = g2p.gid )
group by
  source, evidence_label;

  source         | evidence_label | allele_count
-----------------------+----------------+--------------
cgi                   | B              |            5
cgi                   | C              |           27
cgi                   | D              |          128
civic                 | C              |            4
civic                 | D              |            6
jax                   | C              |            2
jax                   | D              |           21
jax_trials            | D              |            3
molecularmatch_trials | B              |            1
molecularmatch_trials | C              |            2
molecularmatch_trials | D              |            1
oncokb                | D              |           15
(12 rows)

  Time: 64810.378 ms (01:04.810)

```

binary for jsquery
```
# https://wiki.postgresql.org/wiki/Apt
# https://www.ubuntuupdates.org/package/postgresql/trusty-pgdg/main/base/postgresql-11-jsquery


sudo apt-get update
sudo apt-get install -y wget ca-certificates
wget --quiet -O - https://www.postgresql.org/media/keys/ACCC4CF8.asc | sudo apt-key add -
sudo sh -c 'echo "deb http://apt.postgresql.org/pub/repos/apt/ $(lsb_release -cs)-pgdg main" > /etc/apt/sources.list.d/pgdg.list'
sudo apt-get install -y postgresql-10-jsquery

```

Use jsquery @@ syntax to find alleles in the vicinity of chr10:35305276

```
select gid, data->'start' as "start", data->'chromosome' as "chromosome"   from vertex where label = 'Allele' and data @@ 'start ($ > 35305274 AND $ < 35305278) ';

gid                       |  start   | chromosome
-------------------------------------------------+----------+------------
Allele:000010b2803f336a076a45db8cc4ccc34a587680 | 35305277 | "10"
Allele:41f2a16500d995b488a31e4a32f97308ef794149 | 35305277 | "10"
Allele:b0b64df19182d3fe63d9edc3856f97446126f353 | 35305275 | "10"
(3 rows)

Time: 0.968 ms
```


select
  g2p.data->>'source' as "source",
  g2p.data->>'evidence_label' as "evidence_label",
  alleles.project_gid,
  count(distinct alleles.aliquot_gid) as "aliquot_count",
  count(distinct alleles.allele_gid) as "allele_count",
  a.data#>>'{annotations,myvariantinfo,snpeff,ann,0,putative_impact}' as "putative_impact"
from
  alleles  -- edge joins view
    join vertex g2p on (g2p.label = 'G2PAssociation' and alleles.g2passociation_gid = g2p.gid )
    join vertex a on (a.label = 'Allele' and alleles.allele_gid = a.gid )
group by
  source, evidence_label, project_gid, putative_impact;





`Given an Allele, I'd like to find Aliquots, Individuals, Treatments, etc for Genes that are overexpressed`

```
# select * from matrix where key ='ENSG00000000419' and val > (select avg(val)*4 from matrix where key ='ENSG00000000419') ;
            gid             |       key       |       val
----------------------------+-----------------+------------------
 Expression:ccle:ACH-000530 | ENSG00000000419 | 287.675842285156
 Expression:ccle:ACH-000277 | ENSG00000000419 |  218.65934753418
 Expression:ccle:ACH-000050 | ENSG00000000419 | 237.621643066406
(3 rows)

Time: 20.606 ms
```



select * from matrix where  val > (select avg(val)*4 from matrix m where m.key = matrix.key )  limit 1 ;

create table  gene_expression_average  as select key,  avg(val) from matrix group by key  ;
create index gene_expression_average_idx on gene_expression_average (key) ;

select key, val  from matrix where gid = 'Expression:ccle:ACH-000728' and  val > (select a.avg*4 from gene_expression_average a where a.key = matrix.key )   ;
select gid  from vertex where label = 'Allele' and data @@ 'annotations.myvariantinfo.snpeff.ann.#.putative_impact = "HIGH"' limit  1;


It is very expensive to aggregate from callset->allele->gene
```
create or replace view alleles
  as
  select
    p.to      as "project_gid",
    i.to      as "individual_gid",
    i.from    as "biosample_gid",
    al.from   as "allele_gid",
    c.from    as "callset_gid",
    g2p.from  as "g2passociation_gid",
    ph.to     as "phenotype_gid",
    e.to      as "environment_gid",
    t.to      as "treatment_gid",
    g.to      as "gene_gid"
  from
    edge as a
      join edge i on (a.label = 'AliquotFor' and i.from = a.to and i.label = 'BiosampleFor')
        join edge p on (p.label = 'InProject' and p.from = i.to)
        join edge t on (t.label = 'TreatedWith' and t.from = i.to)
          join edge c on (c.label = 'CallsetFor' and a.from = c.to)
            join edge al on (al.label = 'AlleleCall' and c.from = al.to)
              join edge g on (g.label = 'AlleleIn' and al.from = g.from)
              join edge g2p on (g2p.label = 'HasAlleleFeature' and g2p.to = al.from)
                join edge ph on (ph.label = 'HasPhenotype' and g2p.from = ph.from)
                join edge e on (e.label = 'HasEnvironment' and g2p.from = e.from)
                join vertex g2p_v on (g2p_v.label = 'G2PAssociation' and g2p.from = g2p_v.gid)
;
```

add link to expression
```
create or replace view alleles
  as
  select
    p.to      as "project_gid",
    i.to      as "individual_gid",
    i.from    as "biosample_gid",
    c.to      as "aliquot_gid",
    ex.from   as "expression_gid",
    c.from    as "callset_gid",
    al.from   as "allele_gid",
    g2p.from  as "g2passociation_gid",
    ph.to     as "phenotype_gid",
    e.to      as "environment_gid",
    t.to      as "treatment_gid"
  from
    edge as a
      join edge i on (a.label = 'AliquotFor' and i.from = a.to and i.label = 'BiosampleFor')
        join edge p on (p.label = 'InProject' and p.from = i.to)
        join edge t on (t.label = 'TreatedWith' and t.from = i.to)
          join edge c on (c.label = 'CallsetFor' and a.from = c.to)
          join edge ex on (ex.label = 'ExpressionOf' and ex.to = c.to)
            join edge al on (al.label = 'AlleleCall' and c.from = al.to)
              join edge g2p on (g2p.label = 'HasAlleleFeature' and g2p.to = al.from)
                join edge ph on (ph.label = 'HasPhenotype' and g2p.from = ph.from)
                join edge e on (e.label = 'HasEnvironment' and g2p.from = e.from)
                join vertex g2p_v on (g2p_v.label = 'G2PAssociation' and g2p.from = g2p_v.gid)
;
```





select
  alleles.project_gid,
  g2p.data->>'evidence_label' as "evidence_label",
  count(distinct alleles.allele_gid) as "allele_count"
from
  alleles  -- edge joins view
    join vertex g2p on (g2p.label = 'G2PAssociation' and alleles.g2passociation_gid = g2p.gid )
    join vertex a on (a.label = 'Allele' and alleles.allele_gid = a.gid )
group by
  project_gid, evidence_label ;



  select
    alleles.project_gid,
    a.data#>>'{annotations,myvariantinfo,snpeff,ann,0,putative_impact}' as "putative_impact",
    count(distinct alleles.allele_gid) as "allele_count"
  from
    alleles  -- edge joins view
      join vertex a on (a.label = 'Allele' and alleles.allele_gid = a.gid )
  group by
    project_gid, putative_impact;


use jq and psql to load

```
# load a vertex
 head -3  outputs/g2p/Phenotype.Vertex.json  | jq -rc  '"\(.gid)|\(.label)|\(.data)"' | psql -d test -c "copy vertex(gid, label, data) from stdin csv quote '^' delimiter '|' ;"
# load a vertex first by loading all, except a particular json path
zcat  outputs/ccle/Expression.Vertex.json.gz      | jq -rc  '.gid as $gid |  .data.values | to_entries  | map( . + {gid: $gid} )[]  | select(.value != 0 )  |  "\(.gid)|\(.key)|\(.value)"  '     | ~/sql -c "copy matrix(gid, key, val) from stdin csv quote '^' delimiter '|' ;"
# load matrix using gid and data.values transposed ( 1 row per matrix value)
zcat  outputs/ccle/Expression.Vertex.json.gz  | head -1  \
    | jq -rc  '.gid as $gid |  .data.values | to_entries  | map( . + {gid: $gid} )[]  | select(.value != 0 )  |  "\(.gid)|\(.key)|\(.value)"  ' \
    | ~/sql -c "copy matrix(gid, key, val) from stdin csv quote '^' delimiter '|' ;"



cat  outputs/tcga/tcga.Expression.Vertex.json | \
  jq -rc  'del(.data.values) | "\(.gid)|\(.label)|\(.data)"' | \
  ~/sql -c "copy vertex(gid, label, data) from stdin csv quote '^' delimiter '|' ;"  &
cat  outputs/tcga/tcga.Expression.Vertex.json  \
    | jq -rc  '.gid as $gid |  .data.values | to_entries  | map( . + {gid: $gid} )[]  | select(.value != 0 )  |  "\(.gid)|\(.key)|\(.value)"  ' \
    | ~/sql -c "copy matrix(gid, key, val) from stdin csv quote '^' delimiter '|' ;"    &


```
