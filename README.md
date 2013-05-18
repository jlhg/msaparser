# msaparser

Parsing clustal file for sequence variation analysis

## Usage

```python
import msaparser

parser = msaparser.Parse()
parser.parse('protein_sequences.clustal', 'a')
parser.block  # the list of parsed blocks
parser.mutnum  # the number of point mutation
parser.mutprofile  # the list of mutation profiles
parser.mutposlist  # the list of mutation position
parser.get_clustal_html()  # html output of parsed result
```