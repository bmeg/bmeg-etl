

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


select count(distinct v.gid) from vertex as v where v.label = 'Expression'
data.values.ENSG00000227232

select * from vertex as v where v.label = 'Expression' limit 1;


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


CREATE INDEX vertex_ENSG00000227232 ON vertex ((data->'values'->'ENSG00000227232') );



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



CREATE INDEX vertex_data on vertex USING gin (data);
