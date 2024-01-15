import os
import csv
import sys
from code.classes import protein
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
    Main function to execute protein folding based on the specified algorithm and visualize the results.

    The function reads the command line arguments to determine the fold algorithm, dimensions, and
    whether the protein sequence includes type 'C'. It then dynamically imports the specified algorithm,
    applies it to protein sequences from a CSV file, and visualizes the results in 2D or 3D.

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
    The function assumes that the protein sequences are stored in CSV files located in the "data/input/" directory.
    The output CSV files and visualizations will be saved in the "data/output/csv/" and "data/output/plot/"
    directories, respectively.
    """
    if len(sys.argv) != 3 and len(sys.argv) != 4:
        print("Usage: python main.py <fold_algorithm> <dimensions> [C]")
        sys.exit(1)

    fold_algorithm: str = sys.argv[1].lower()

    if fold_algorithm not in ALGORITHM_FILES:
        raise ValueError("Invalid fold type.")

    # Import the selected folding algorithm dynamically
    fold_module = __import__(
        f"code.algorithms.{fold_algorithm}", fromlist=[fold_algorithm])
    fold_function = getattr(fold_module, f"{fold_algorithm}_fold")

    if len(sys.argv) == 3:
        filename: str = "data/input/sequences_H_P.csv"

    if len(sys.argv) == 4:
        if sys.argv[3] != "C":
            print("Usage: python main.py <fold_algorithm> <dimensions> [C]")
            sys.exit(2)

        filename: str = "data/input/sequences_H_P_C.csv"

    with open(filename, "r") as file:
        reader = csv.reader(file)
        line_number: int = 0

        for row in reader:
            if not row:
                break

            sequence: str = row[0]
            test_protein: protein.Protein = protein.Protein(sequence)

            # Call the selected folding algorithm
            fold_function(test_protein)

            test_protein.create_csv(fold_algorithm, line_number)
            if sys.argv[2] == "2":
                visualization_2D.plot_2d(
                    test_protein, ("red", "blue", "green"), line_number, fold_algorithm)
            elif sys.argv[2] == "3":
                visualization_3D.plot_3d(
                    test_protein, ("red", "blue", "green"), line_number, fold_algorithm)
            else:
                print("Please enter dimension as 2 or 3")
                sys.exit(3)

            line_number += 1
            # break is for test (1x plotten)
            # break


if __name__ == "__main__":
    main()
