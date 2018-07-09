
class: Workflow
cwlVersion: v1.0

$namespaces:
  bmeg: http://bmeg.io

bmeg:split: source/pubmed.list

inputs:
  url:
    type: string
    bmeg:input: url

requirements:
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement

steps:
  download:
    run: curl/curl.cwl
    in:
      URL: url
      NAME:
        valueFrom: $(inputs.URL.split("/").pop())
    out:
      - OUTPUT

outputs:
  OUT:
    type: File
    outputSource: download/OUTPUT
    bmeg:key: source/pubmed/{name}
