import csv
import sys
from code.classes import protein
from code.visualization import visualization_2D
from code.visualization import visualization_3D


def main():
    showGrid = True if len(sys.argv) == 2 and sys.argv[1] == "grid" else False

    with open("data/input/sequences_H_P.csv", "r") as file:
        reader = csv.reader(file)
        line_number = 0

        for row in reader:
            if not row:
                break

            sequence = row[0]
            test_protein = protein.Protein(sequence)
            test_protein.create_csv(line_number)
            visualization_3D.plot_3d(test_protein, ("red", "blue", "green"))
            visualization_2D.plot_2d(test_protein, ("red", "blue", "green"))

            if showGrid:
                print(test_protein.get_grid_2D())

            line_number += 1
            # break is for test (1x plotten)
            break


if __name__ == "__main__":
    main()
