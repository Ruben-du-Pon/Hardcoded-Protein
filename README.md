# Hardcoded Protein

## Usage

First install requirements using

```bash
pip install -r requirements.txt
```

Then run the program by running

```bash
python main.py <fold_algorithm> <dimensions> <iterations> [png/svg] [C/c]
```

where `<fold_algorithm>` is the filename (without .py) of one of the algorithms in the code/algorithms folder, `<dimensions>` is either 2 or 3 for 2D or 3D folding, respectively, `<iterations>` is the number of times the algorithm should be run, the `[png/svg]` is to specify the output format and the optional argument `[C]` adds the Cysteine aminoacid proteins.


Below follows a table on the estimate time to run a certain algorithm in 2D or 3D based on the length of the protein itself.

## 2D

```
|        Algorithm         |  Estimated time, length ~10   |  Estimated time, length ~10-20   |      Avg length ~10      |      Avg length ~10-20      |
| -------------------------| ----------------------------- | -------------------------------- | ------------------------ | --------------------------- |
| Baseline                 |               ~               |              ~                   |            ~             |             ~               |
| Simulated Annealing      |               ~               |              ~                   |            ~             |             ~               |
| BFS with Random sampling |               ~               |              ~                   |            ~             |             ~               |
| FRESS                    |               ~               |              ~                   |            ~             |             ~               |
```


## 3D

```
|        Algorithm         |  Estimated time, length ~10   |  Estimated time, length ~10-20   |      Avg length ~10      |      Avg length ~10-20      |
| -------------------------| ----------------------------- | -------------------------------- | ------------------------ | --------------------------- |
| Baseline                 |               ~               |              ~                   |            ~             |             ~               |
| Simulated Annealing      |               ~               |              ~                   |            ~             |             ~               |
| BFS with Random sampling |               ~               |              ~                   |            ~             |             ~               |
| FRESS                    |               ~               |              ~                   |            ~             |             ~               |

```
