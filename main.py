import os
import csv
import sys
from code.classes.protein import Protein
from code.visualization import visualization_2D
from code.visualization import visualization_3D

# Dynamically get filenames without extensions from the algorithms directory
ALGORITHM_FILES = [
    os.path.splitext(file)[0]
    for file in os.listdir("code/algorithms")
    if file.endswith(".py") and file != "__init__.py"
]


def main() -> None:
    """
    Main function to execute protein folding based on the specified algorithm
    and visualize the results.

    The function reads the command line arguments to determine the fold
    algorithm, dimensions, whether the protein sequence includes type 'C', and
    the output file format. It then dynamically imports the specified algorithm,
    applies it to protein sequences from a CSV file, and visualizes the results in
    2D or 3D.

    Parameters
    ----------
    None

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If an invalid fold type is provided.

    Notes
    -----
    The function assumes that the protein sequences are stored in CSV files
    located in the "data/input/" directory.
    The output files will be saved in the "data/output/csv/" and "data/output/plot/"
    directories for CSV and visualization files, respectively.
    """
    filename: str = "data/input/sequences_H_P.csv"

    if len(sys.argv) != 4 and len(sys.argv) != 5 and len(sys.argv) != 6:
        print(
            "Usage: python main.py <fold_algorithm> <dimensions> <iterations> [png/svg] [C/c]")
        sys.exit(1)

    fold_algorithm: str = sys.argv[1].lower()

    if sys.argv[2] not in ("2", "3"):
        print("Please enter dimension as 2 or 3")
        sys.exit(2)

    dimensions: int = int(sys.argv[2])

    if sys.argv[3].isdigit() is False:
        print("Please enter iterations as an integer")
        sys.exit(3)

    iterations: int = int(sys.argv[3])

    # Check if the specified fold algorithm is valid
    if fold_algorithm not in ALGORITHM_FILES:
        raise ValueError("Invalid fold type.")

    # Import the selected folding algorithm dynamically
    fold_module = __import__(f"code.algorithms.{fold_algorithm}", fromlist=[
                             f"{fold_algorithm}Fold"])
    fold_class = getattr(fold_module, f"{fold_algorithm.capitalize()}Fold")

    # Read the protein sequences from the CSV file
    if len(sys.argv) == 4:
        filename: str = "data/input/sequences_H_P.csv"

    if len(sys.argv) == 5:
        if sys.argv[4].lower() not in ("png", "svg"):
            print("Please enter file format as png or svg")
            sys.exit(5)

    if len(sys.argv) == 6:
        if sys.argv[5].lower() == "c":
            filename: str = "data/input/sequences_H_P_C.csv"

    with open(filename, "r") as file:
        reader = csv.reader(file)
        line_number: int = 0

        # Read the file line by line until an empty line is reached
        for row in reader:
            if not row:
                break

            # Create Protein objects for each sequence
            sequence: str = row[0]
            test_protein: Protein = Protein(sequence)

            # Initialize the folding algorithm class
            fold_instance = fold_class(test_protein, dimensions, iterations)

            # Call the run method of the folding algorithm
            fold_instance.run()

            # Set the output filenames
            if dimensions == 2:
                filename = f"data/output/csv/{fold_algorithm}_{line_number}_2D.csv"
                plotname = f"data/output/plot/{fold_algorithm}_{line_number}_2D.{sys.argv[4].lower() if len(sys.argv) == 5 else 'svg'}"
            else:
                filename = f"data/output/csv/{fold_algorithm}_{line_number}_3D.csv"
                plotname = f"data/output/plot/{fold_algorithm}_{line_number}_3D.{sys.argv[4].lower() if len(sys.argv) == 5 else 'svg'}"

            # Write the results to a CSV file and the visualization to a file
            test_protein.create_csv(filename)
            if dimensions == 2:
                visualization_2D.plot_2d(
                    test_protein, ("red", "blue", "green"), plotname, sys.argv[4].lower() if len(sys.argv) == 5 else 'svg')
            else:
                visualization_3D.plot_3d(
                    test_protein, ("red", "blue", "green"), plotname, sys.argv[4].lower() if len(sys.argv) == 5 else 'svg')

            line_number += 1


if __name__ == "__main__":
    main()
