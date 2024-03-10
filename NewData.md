
# Adding New Data to BMEG

## Processing and data dependancies

Add build files in BMEG are organized using [Lathe](https://github.com/bmeg/lathe). All `.plan` 
files in the import module will be scanned to determine which command lines need to be run in order 
to build dependencies. 

### Example Plan file
Lathe files are written in Javascript. The program is executed to declare workflows and processes. 
Each process is declared with the following properties:

| commandLine | string | Command line to be invoked |
| shell  |  string | Similar to `CommandLine` except executed using BASH to evaluate in scripting commands, such as pipe commands `>` | 
| inputs | map[string]string | Key/Value mapping of named inputs to their paths |
| outputs | map[string]string |
| memMB  | uint | Number of MB of RAM required for the job |
| nppus  | uint | Number of CPUs to be reseved for the job |

The `commandLine` and `shell` commands are evaluated using the [Mustache](https://mustache.github.io/) templating language, 
with the `inputs` and `outputs` maps being mad avalible. 


Example script:
```javascript
prep = lathe.Workflow("prep")

prep.Add(lathe.Process({ 
    commandLine: "curl https://ftp.expasy.org/databases/cellosaurus/cellosaurus.obo -o {{outputs.cellObo}}",
    outputs: {
       cellObo: "../../source/cellosaurus/cellosaurus.obo"
    }
}))
```



## Sifter files
[Sifter](https://bmeg.github.io/sifter/) is a stream based processing program and can be used to transform data.
Sifter uses YAML based declaration files, that allow the management system to easily scan the flow of data 
and determine where different types of data being stored. 


## Object schema validation and storage
To identify formatted data that will be incoprated into the evidence graph, Sifter looks for pipelines that run 
object validation and then emit objects to files. An example of this two-part stanza in the pipeline would be:

```yaml
    - objectValidate:
        title: Publication
        schema: "{{config.schema}}"
    - emit:
        name: Publication
```

Now the management system will know that correctly formatted `Publication` objects will be found in the file created by the `emit` step.

### Sifter output files

Emitted sifter output files are written into the following path format:
```
<outdir>/<plan name>.<pipeline name>.<emit name>.json.gz
```

So the following plan file would place output files in `../../output/pubmed/pubmed.transform.Publication.json.gz`

```yaml
class: sifter
name: pubmed
outdir: ../../output/pubmed

config:
  schema: ../../schema

...

pipelines:
  transform:
    ...
    - objectValidate:
        title: Publication
        schema: "{{config.schema}}"

    - emit:
        name: Publication

```

### Installing required tools

Install sifter 

```
go install github.com/bmeg/sifter
```

Install lathe
```
go install github.com/bmeg/lathe
```


To invoke a Lathe file directly
```
lathe run ./path/to/build.plan
```

To invoke a Sifter file directly
```
sifter run ./path/to/sifter.yaml
```


To generate graph build steps

```
sifter graph-plan -o ./graph -C ./graph-build path/to/sifter.yaml
```

This command generates the scripts to build the graph elements, 
but does not actually build the graph. To build the graph, the sifter
command needs to be run on the plan file that was generated. 
The `-o` flag determines where graph elements (node and edge files) will
be placed when the graph build, the `-C` file changes the directory that 
generated sifter file will be placed. 

So the command
```
sifter graph-plan -o graph -C graph-build transform/pubmed/transform.yaml
```
will generate the file

```yaml
name: pubmed-graph
class: sifter

outdir: ../graph

config:
    transform: ../output/pubmed/pubmed.transform.publication.json.gz
    transformSchema: ../schema
inputs:
    transform:
        jsonLoad:
            input: "{{config.transform}}"
pipelines:
    transform-graph:
        - from: transform
        - graphBuild:
            schema: "{{config.transformSchema}}"
            title: Publication
```

When the command `sifter run graph-build/pubmed.yaml` is run it will generate the files
```
graph/pubmed-graph.transform-graph.vertex.json.gz
graph/pubmed-graph.transform-graph.edge.json.gz
```


### Loading graph data into a test GRIP instance

Install GRIP
```
go get github.com/bmeg/grip
```

Turn grip of with Pebble driver
```
grip server -d pebble
```

Create a graph
```
grip create test-db
```

Load database
```
grip load test-db --dir test-graph/
```

Query the graph in python:
```
import gripql
conn = gripql.Connection("http://localhost:8201")
G=conn.graph("test-db")
G.query().V().count().execute()
```