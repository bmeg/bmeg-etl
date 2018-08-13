### Directory layout

The source, transform, outputs, and tests trees should mirror each other:

```
├── outputs
│   ├── ccle
│   ├── gdsc
│   ├── gtex
├── source
│   ├── ccle
│   ├── gdsc
│   ├── gtex
├── tests
│   ├── ccle
│   ├── gdsc
│   ├── gtex
├── transform
│   ├── ccle
│   ├── gdsc
│   ├── gtex
```

Transformers should read from `source` and write to `outputs`. Directory names should reflect the project or data source.

### Transformers should be specific

All transformers should run without CLI arguments, e.g. `python transform/mc3/variants.py`.

Transformers should be specific to one data source, e.g. `transform/mc3` or `transform/ccle`.

Code should expect to be run from the root of the repo. Using the directory layout defined above, transformers should
open files such as `gzip.open("source/mc3/variants.maf.gz")`. Emitters and caching utils (requests) should handle
output paths automatically.

This allows us to collaborate more easily, quickly, consistently, and correctly. We no longer need to understand the correct
CLI arguments to use when working on someone else's transformer. It also documents and commits/versions the arguments
used by including it in the python code.

### End-to-end tests

End-to-end tests should be included in the `tests` directory. These tests should make AQL queries against Arachne
and verify the results.

### HTTP requests

All HTTP requests should use `bmeg.requests.Client`, if possible. This client is preconfigured with retries, timeouts,
caching, and more.

### Python

Python 3.7+ only

### Gzip files

This pattern is extremely common:
```
if path.endswith(".gz"):
  gzip.open(path)
else:
  open(path)
```

Use this instead: `bmeg.ioutils.reader(path)`

### GID constructors

Where possible, GIDs should be a simple join of model fields. For example:

```
class Foo:
  bar: str

  def make_gid(self, bar):
    return "Foo:" + bar.
```

### HG19 genome only 

All data should be based on HG19 genome reference.

### 0-based, half-open genomic intervals

All genomic intervals should be 0-based, half-open.

### Don't over generalize

Be careful with creating reusable transformer code. Almost all inputs have some slight, subtle variation.
Reusable code for these cases can easily end up adding a ton of options to allow handling of all these cases.
In the end, such code becomes complex and fragile.

It's better to have simple code and a little bit of duplication.

### Model naming

Vertex model names are `UpperCamelCase`. Edge model names are `lowerCamelCase`.  
All field names are `lower_snake_case`.

### Edge direction.

Edges should point from child to parent, if possible.  
For example, `Exon -> Transcript -> Gene` not `Gene -> Transcript -> Exon`.

### Model field types

Transformer output (i.e. JSON) should match the defined model fields, if possible.  
If a field might be None/null, that should be defined on the model, e.g.

```python
from typing import Union

class DrugResponse:
  # Some samples have a "NA" value for drug response
  value: Union[None, float]
```

### Code style: remove dead/commented code

Dead code is confusing and a maintenance burder. If code is commented out, it's not necessary and should be removed.
Same for dead code (functions never called, code after `return` statements, etc).

### Code style: no shebags or magic comments

```
#!/usr/bin/env python
# -*- coding: utf-8 -*-
```

Don't add these.

Everyone has some different version of these. None of them are necessary. Some of them are wrong.

### Code style: don't use the _private_var pattern

```python
def _private_func(arg): pass
```

It's very easy to create inconsistent code with this pattern, which is confusing. We're not writing library code,
so there's no need for private variables/functions.

### Code style: don't bypass the linter

```python
raise Exception('something failed and I need to tell you about it {}'.format(var1))  # noqa
```

Don't use the `# noqa` escape hatch. There's no point to code linting if it can be skipped whenever we feel like it.
If our linting rules are annoying, we should change them, otherwise we should stick to them.
