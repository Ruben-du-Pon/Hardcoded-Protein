import cProfile
import sys
from importlib import import_module
from codefiles.classes.protein import Protein


def main() -> None:
    """
    Entry point of the profiling script.
    Parses command line arguments, imports the specified algorithm module,
    and runs the algorithm with the given number of iterations.

    Args:
        None

    Returns:
        None
    """
    argc = len(sys.argv)
    argv = sys.argv
    if argc != 3:
        print("Usage: python3 profiling.py <algorithm> <iterations>")
        exit(1)

    if not argv[2].isdigit():
        raise ValueError("Iterations must be an integer.")

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
