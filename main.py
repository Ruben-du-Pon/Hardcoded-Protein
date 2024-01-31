import os
import sys
import subprocess

# Dynamically get filenames without extensions from the algorithms directory
ALGORITHM_FILES = [
    os.path.splitext(file)[0]
    for file in os.listdir("codefiles/algorithms")
    if file.endswith(".py") and file != "__init__.py"
]


def main() -> None:
    """
    Main function that reads protein sequences from a file and generates protein folds using the specified algorithm.

    Usage: python main.py <fold_algorithm> <dimensions> <iterations> [C/c]

    Parameters:
    - fold_algorithm (str): The algorithm to use for protein folding.
    - dimensions (int): The number of dimensions for the protein fold (2 or 3).
    - iterations (int): The number of iterations for the protein folding process.
    - C/c (optional): Specify 'C' or 'c' to use a different input file.

    Raises:
    - ValueError: If the specified fold algorithm is invalid.
    """

    # Check if the correct number of command line arguments is given
    if len(sys.argv) != 3 and len(sys.argv) != 4:
        print(
            "Usage: python main.py <fold_algorithm> <dimensions> [C/c]")
        sys.exit(1)

    fold_algorithm: str = sys.argv[1].lower()

    # Check if the specified dimensions are valid
    if sys.argv[2] not in ("2", "3"):
        print("Please enter dimension as 2 or 3")
        sys.exit(2)

    dimensions: int = int(sys.argv[2])

    # Check if the specified fold algorithm is valid
    if fold_algorithm not in ALGORITHM_FILES:
        raise ValueError("Invalid fold type.")

    # Get a copy of the current environment variables
    env = os.environ.copy()
    # Add the parent directory of codefiles to the Python path in the new environment
    env["PYTHONPATH"] = os.path.abspath(os.path.join(
        os.path.dirname(__file__), '..')) + ":" + env.get("PYTHONPATH", "")

    # Call the generate_<algorithm>.py script with the correct parameters
    subprocess.run(["python3",
                    f"codefiles/fold_generation/generate_{fold_algorithm}.py",
                    str(dimensions), "t" if len(sys.argv) == 4 else "f"],
                   env=env)  # Pass the modified environment to the new process


if __name__ == "__main__":
    main()
