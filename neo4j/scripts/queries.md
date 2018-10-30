

```
cypher-shell -u neo4j -p password

// create a test node from a param

:param foo => { a: "A" }
CREATE(a:Foo $foo) ;

match (m) return m.a ;
+-----+
| m.a |
+-----+
| "A" |
+-----+

:param props => { props : { name : "Andy", position : "Developer" } }
:param props => { name : "Andy", position : "Developer" }


// note file: urls need to point at import dir, defaults to  /var/lib/neo4j/import/
zcat  /src/outputs/ccle/Aliquot.Vertex.json.gz | head -3  > /var/lib/neo4j/import/Aliquot.Vertex.json


// import, just show data
WITH "file:///Aliquot.Vertex.json" AS url
CALL apoc.load.json(url)  YIELD value
RETURN value ;
+------------------------------------------------------------------------------------------------------------+
| value                                                                                                      |
+------------------------------------------------------------------------------------------------------------+
| {_id: "Aliquot:ACH-000557", gid: "Aliquot:ACH-000557", label: "Aliquot", data: {aliquot_id: "ACH-000557"}} |
| {_id: "Aliquot:ACH-001000", gid: "Aliquot:ACH-001000", label: "Aliquot", data: {aliquot_id: "ACH-001000"}} |
| {_id: "Aliquot:ACH-000198", gid: "Aliquot:ACH-000198", label: "Aliquot", data: {aliquot_id: "ACH-000198"}} |
+------------------------------------------------------------------------------------------------------------+


// import, set property(ies)
// this works because all members of Aliqout's 'data' are primitives
WITH "file:///Aliquot.Vertex.json" AS url
  CALL apoc.load.json(url) YIELD value
      CREATE (t:Test)
        SET t = value.data,
            t.gid = value.gid
;

match (m:Test) return ( m ) ;
+---------------------------------------------------------------+
| m                                                             |
+---------------------------------------------------------------+
| (:Test {aliquot_id: "ACH-000557", gid: "Aliquot:ACH-000557"}) |
| (:Test {aliquot_id: "ACH-001000", gid: "Aliquot:ACH-001000"}) |
| (:Test {aliquot_id: "ACH-000198", gid: "Aliquot:ACH-000198"}) |
+---------------------------------------------------------------+



// import, set property(ies)
// this fails because  members of Allele's 'data' are primitives & maps
WITH "file:///Allele.Vertex.json" AS url
  CALL apoc.load.json(url) YIELD value
      CREATE (t:Test)
        SET t = value.data,
            t.gid = value.gid
;

// solution is to filter the data map removing maps and lists
// load a file of json
WITH "file:///Allele.Vertex.json" AS url
  // yield each line
  CALL apoc.load.json(url) YIELD value
    // create a vertex with dynamic label
    // set its values, filtering out lists and maps (only primitives allowed)    
    CALL apoc.create.node(
      [value.label],
      apoc.map.clean(value.data, filter(k IN keys(value.data) WHERE apoc.meta.type(value.data[k]) in ['MAP', 'LIST']), [])
    ) YIELD node
    // include the vertex gid
    SET node.gid = value.gid
    // return the count of nodes created
    RETURN count(node) as nodes_created
;

match(a:Allele) return (a) ;
+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| a                                                                                                                                                                               |
+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| (:Allele {genome: "GRCh37", gid: "Allele:bef7d80d2d7a1f63a4be8dc7d43b8c867f1c29e2", alternate_bases: "T", chromosome: "1", start: 1277461, end: 1277461, reference_bases: "C"}) |
| (:Allele {genome: "GRCh37", gid: "Allele:f14ebbb5eca5b03314d6d5773588c256bd6bcdea", alternate_bases: "A", chromosome: "1", start: 2144416, end: 2144416, reference_bases: "G"}) |
| (:Allele {genome: "GRCh37", gid: "Allele:4cb95b67ef30051e36ab4f69aa7d728262eab76c", alternate_bases: "C", chromosome: "1", start: 2435359, end: 2435359, reference_bases: "A"}) |
+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+


// solution is to filter the data map removing maps and lists

// load a file of json
WITH "file:///outputs/ccle/Allele.Vertex.json.gz" AS url
  // yield each line
  CALL apoc.load.json(url) YIELD value
    WITH value limit 1
      // create a vertex with dynamic label
      // set its values, filtering out lists and maps (only primitives allowed)    
      CALL apoc.create.node(
        [value.label],
        apoc.map.clean(value.data, filter(k IN keys(value.data) WHERE apoc.meta.type(value.data[k]) in ['MAP', 'LIST']), [])
      ) YIELD node
      // include the vertex gid
      SET node.gid = value.gid
      // return the count of nodes created
      RETURN count(node) as nodes_created
;


// zipped files work too!
WITH "file:///outputs/ccle/Aliquot.Vertex.json.gz" AS url
  // yield each line
  CALL apoc.load.json(url) YIELD value
    WITH value LIMIT 1
      // create a vertex with dynamic label
      // set its values, filtering out lists and maps (only primitives allowed)    
      CALL apoc.create.node(
        [value.label],
        apoc.map.clean(value.data, filter(k IN keys(value.data) WHERE apoc.meta.type(value.data[k]) in ['MAP', 'LIST']), [])
      ) YIELD node
      // include the vertex gid
      SET node.gid = value.gid
      // return the count of nodes created
      RETURN count(node) as nodes_created
;


// edges


// print first edge
WITH "file:///outputs/ccle/AliquotFor.Edge.json.gz" AS url
  // yield each line
  CALL apoc.load.json(url) YIELD value
    RETURN  value LIMIT 1


WITH "file:///outputs/ccle/AliquotFor.Edge.json.gz" AS url
  // yield each line
  CALL apoc.load.json(url) YIELD value
    WITH value LIMIT 1
      MATCH(f {gid: value.from})
      MATCH(t {gid: value.to})
      // create an vertex with dynamic label
      // set its values, filtering out lists and maps (only primitives allowed)    
      CALL apoc.create.relationship(
        f,
        value.label,
        apoc.map.clean(value.data, filter(k IN keys(value.data) WHERE apoc.meta.type(value.data[k]) in ['MAP', 'LIST']), []),
        t
      ) YIELD rel
      // skip the edge gid, no use case
      // SET rel.gid = value.gid
      // return the count of relationships created
      RETURN count(rel) as relationships_created
;







```

cd /var/lib/neo4j/plugins
wget https://github.com/neo4j-contrib/neo4j-apoc-procedures/releases/download/3.4.0.3/apoc-3.4.0.3-all.jar



CALL apoc.config.list() yield key, value
where key =~ '.*apoc.*'
return key, value ;


CALL apoc.config.list() yield key, value
return key, value ;


CALL apoc.config.list() yield key, value
where key =~ '.*dbms.*'
return key, value ;
