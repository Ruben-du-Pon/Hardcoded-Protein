import cProfile
import sys
from importlib import import_module
from codefiles.classes.protein import Protein


def main():
    argc = len(sys.argv)
    argv = sys.argv
    if argc != 3:
        print("Usage: python3 profiling.py <algorithm> <iterations>")
        exit(1)

    algorithm = str(argv[1]).lower()
    module_name = f"codefiles.algorithms.{algorithm}"
    module = import_module(module_name)

    iterations = int(argv[2])

    class_name = algorithm.capitalize() + "Fold"
    AlgorithmFold = getattr(module, class_name)
    protein = Protein("HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH")

    def run_algorithm():
        # Create an instance of the class
        instance = AlgorithmFold(protein, 2, iterations)
        instance.run()

    output = f"data/output/profiles/{algorithm}.prof"
    cProfile.runctx('run_algorithm()', globals(), locals(), output)


if __name__ == "__main__":
    main()