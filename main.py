import os
import sys
import subprocess

from codefiles.classes.protein import Protein

# Dynamically get filenames without extensions from the algorithms directory
ALGORITHM_FILES = [
    os.path.splitext(file)[0]
    for file in os.listdir("codefiles/algorithms")
    if file.endswith(".py") and file != "__init__.py"
]


def main() -> None:
    """
    Main function that reads protein sequences from a file and generates protein folds using the specified algorithm.

    Usage: python main.py <fold_algorithm> <dimensions> <iterations> <y/n> ['generate']

    Parameters:
    - fold_algorithm (str): The algorithm to use for protein folding.
    - dimensions (int): The number of dimensions for the protein fold (2 or 3).
    - iterations (int): The number of iterations for the protein folding process.
    - C/c (optional): Specify 'C' or 'c' to use a different input file.

    Raises:
    - ValueError: If the specified fold algorithm is invalid.
    """  # noqa

    # Check if the correct number of command line arguments is given
    if len(sys.argv) != 5 and len(sys.argv) != 6:
        print(
            "Usage: python main.py <fold_algorithm> <dimensions> <iterations> <y/n> ['generate']")  # noqa
        sys.exit(1)

    if len(sys.argv) == 6 and sys.argv[5].lower() != "generate":
        print(
            "Usage: python main.py <fold_algorithm> <dimensions> <iterations> <y/n> ['generate']")  # noqa
        sys.exit(1)

    fold_algorithm: str = sys.argv[1].lower()

    # Check if the specified dimensions are valid
    if sys.argv[2] not in ("2", "3"):
        raise ValueError("Invalid dimensions.")

    dimensions: int = int(sys.argv[2])

    # Check if the specified fold algorithm is valid
    if fold_algorithm not in ALGORITHM_FILES:
        raise ValueError("Invalid fold type.")

    # Check if the specified number of iterations is valid
    if not sys.argv[3].isdigit():
        raise ValueError("Invalid number of iterations.")

    iterations = int(sys.argv[3])

    # Check if specification for including the C acid is valid
    if sys.argv[4].lower() not in ("y", "n"):
        raise ValueError("Specify 'y' or 'n' for including the C acid.")

    if len(sys.argv) == 6 and sys.argv[5].lower() == "generate":
        # Get a copy of the current environment variables
        env = os.environ.copy()
        # Add the parent directory of codefiles to the Python path in the new
        # environment
        env["PYTHONPATH"] = os.path.abspath(os.path.join(
            os.path.dirname(__file__), '..')) + ":" + env.get("PYTHONPATH", "")

        # Call the generate_<algorithm>.py script with the correct parameters
        subprocess.run(["python3",
                        f"codefiles/fold_generation/generate_{fold_algorithm}.py",  # noqa
                        str(dimensions),
                        "t" if sys.argv[4].lower == "y" else "f"],
                       env=env)

    elif len(sys.argv) == 5:

        # Determine the file to read the protein sequence from
        if sys.argv[4].lower() == "n":
            sequence_file = "data/input/sequences_H_P"
        else:
            sequence_file = "data/input/sequences_H_P_C"

        # Read the protein sequence from the file
        with open(sequence_file) as f:
            protein_sequence = f.read().strip()

        # Create a Protein object
        protein = Protein(protein_sequence)

        # Dynamically import the specified fold algorithm
        fold_algorithm_module = __import__(
            f"codefiles.algorithms.{fold_algorithm}",
            fromlist=[fold_algorithm])

        # Construct the class name
        class_name = fold_algorithm.capitalize() + "Fold"

        # Get the specified fold algorithm class
        fold_algorithm_class = getattr(
            fold_algorithm_module, class_name)

        # Create an instance of the specified fold algorithm class
        fold_algorithm_instance = fold_algorithm_class(protein,
                                                       dimensions,
                                                       iterations)

        # Call the run method of the specified fold algorithm instance
        fold_algorithm_instance.run()


if __name__ == "__main__":
    main()
