

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
;

allele | callset | aliquot | biosample | individual | project | g2p_association | phenotype | environment | treatment
--------+---------+---------+-----------+------------+---------+-----------------+-----------+-------------+-----------
   174 |    1878 |    1878 |      1878 |        908 |      19 |             102 |       163 |         115 |       150
(1 row)

Time: 24250.170 ms (00:24.250)

```
