

# Common arguments and logging setup


```

from bmeg.util.cli import default_argument_parser
from bmeg.util.logging import default_logging

...

parser = default_argument_parser()
parser.add_argument('--maf_file', type=str,
                    help='Path to the maf you want to import')
# We don't need the first argument, which is the program name
options = parser.parse_args(sys.argv[1:])
default_logging(options.loglevel)

...

```
