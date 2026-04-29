# GBIFOpenTree



Generate files for a Darwin Core Archive file of specimen data
from OpenTree's phylesystem data store

usage:

```python phylesystem_dwca_translator.py --study_id ot_1003```

Generates Nexus fromat tree files for all trees in the study.
Generates an occurence.csv file mapping tip labels to indentifiers and any other available information

Uses python-opentree and dendropy

To setup in a virtual env, run:
```
virtualenv -p python3 venv-gbif-opentree
source venv-gbif-opentree/bin/activate
pip install -r requirements.txt 
```