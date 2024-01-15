import os
import csv
import sys
from code.classes import protein
from code.visualization import visualization_2D
from code.visualization import visualization_3D

# Dynamically get filenames without extensions from the algorithms directory
algorithm_files = [
    os.path.splitext(file)[0]
    for file in os.listdir("code/algorithms")
    if file.endswith(".py") and file != "__init__.py"
]


def main() -> None:
    if len(sys.argv) != 3 or len(sys.argv) != 4:
        print("Usage: python main.py <fold_algorithm> <dimensions> [C]")
        sys.exit(1)

    fold_type = sys.argv[1].lower()

    if fold_type not in algorithm_files:
        raise ValueError("Invalid fold type.")

    # Import the selected folding algorithm dynamically
    fold_module = __import__(
        f"code.algorithms.{fold_type}", fromlist=[fold_type])
    fold_function = getattr(fold_module, f"{fold_type}_fold")

    if len(sys.argv) == 3:
        filename = "data/input/sequences_H_P.csv"

    if len(sys.argv) == 4:
        if sys.argv[3] != "C":
            print("Usage: python main.py <fold_algorithm> <dimensions> [C]")
            sys.exit(2)

        filename = "data/input/sequences_H_P_C.csv"

    with open(filename, "r") as file:
        reader = csv.reader(file)
        line_number = 0

        for row in reader:
            if not row:
                break

            sequence = row[0]
            test_protein = protein.Protein(sequence)

            # Call the selected folding algorithm
            fold_function(test_protein)

            test_protein.create_csv(line_number)
            if sys.argv[2] == "2":
                visualization_3D.plot_2d(
                    test_protein, ("red", "blue", "green"), line_number)
            elif sys.argv[2] == "3":
                visualization_2D.plot_3d(
                    test_protein, ("red", "blue", "green"), line_number)
            else:
                print("Please enter dimension as 2 or 3")
                sys.exit(3)

            line_number += 1
            # break is for test (1x plotten)
            # break


if __name__ == "__main__":
    main()
