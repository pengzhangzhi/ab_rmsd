# AB RMSD
![cover](assets/cover.png)

## Table of Contents

- [About](#about)
- [Getting Started](#getting_started)
- [Usage](#usage)

## About <a name = "about"></a>

Calculate the RMSD between two antibody structure (including nanobody and antibody).


## Getting Started <a name = "getting_started"></a>


First, install abnumber using conda.
```bash
conda install -c bioconda abnumber
```
Then, install ab_rmsd from github.
```bash
pip install git+https://github.com/pengzhangzhi/ab_rmsd.git
```
#### [Optional] DockQ
If you want to calculate the DockQ between complex structures, run the following command.
```bash
cd DockQ
make
```

## Usage <a name = "usage"></a>

Calculate the RMSD between predicted and native nanobody structure.
```python
from ab_rmsd import calc_ab_rmsd

# two example pdb files are provided in the `exampl`e folder.
native = "example/7d6y_1_B.pdb"
pred = "example/pred_7d6y_1_B.pdb"
rmsd = calc_ab_rmsd(native,pred)
print(rmsd)

"""
output:
{'CDRH1': tensor(0.3136), 'CDRH2': tensor(0.0898), 'CDRH3': tensor(4.3704), 'fv-H': tensor(0.6426)}
"""
``` 

# TODOs
- add `DockQ` as an evaluation for heavy and light chain complex structure.

# Credits

Part of the code is adapted from [shitong's Diffab](https://github.com/luost26/diffab).
