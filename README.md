# Hardcoded Protein

## Usage

First install requirements using

```bash
pip install -r requirements.txt
```

Then run the program by running

```bash
python main.py <fold_algorithm> <dimensions> <iterations> [C]
```

where `<fold_algorithm>` is the filename (without .py) of one of the algorithms in the code/algorithms folder, `<dimensions>` is either 2 or 3 for 2D or 3D folding, respectively, `<iterations>` is the number of times the algorithm should be run and the optional argument `[C]` adds the
Cysteine aminoacid proteins.
