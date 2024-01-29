import csv
import sys
from codefiles.algorithms.annealing import AnnealingFold
from codefiles.classes.protein import Protein


def generate_data(dimensions: int, C: bool) -> None:
    """
    Generate a baseline for a protein folding problem.

    Parameters
    ----------
    dimensions : int
        The number of dimensions in which the protein folding will be simulated.
    C : bool
        A boolean indicating whether the protein sequence includes cysteine (C) or not.
    """  # noqa
    # Set the input and output filenames
    if C:
        filename: str = "data/input/sequences_H_P_C.csv"
        outputfile: str = f"data/output/annealing_data/{dimensions}D_C.csv"
    else:
        filename: str = "data/input/sequences_H_P.csv"
        outputfile: str = f"data/output/annealing_data/{dimensions}D.csv"

    # Empty the output file
    with open(outputfile, "w") as file:
        pass

    with open(filename, "r") as file:
        reader = csv.reader(file)
        line_number: int = 0

        # Read the file line by line until an empty line is reached
        for row in reader:
            if not row:
                break

            # Create Protein objects for each sequence
            sequence: str = row[0]

            scores = []

            test_protein: Protein = Protein(sequence)
            test = AnnealingFold(
                test_protein, dimensions, 2500, scores, outputfile)
            test_protein = test.run()
            scores = test.get_scores()

            # Write the average score to the output file
            with open(outputfile, "a") as file:
                writer = csv.writer(file)
                writer.writerow([])
                writer.writerow(
                    ["Average:", f"{sum(scores) / len(scores)}"])
                writer.writerow([])
                writer.writerow([])

            line_number += 1


def main() -> None:
    """
    The main function that checks command line arguments and calls the
    generate_data function.
    """
    if len(sys.argv) != 3:
        print("Usage: python generate_hillclimber_data.py <dimensions> <t/f>")
        sys.exit(1)

    dimensions: int = int(sys.argv[1])
    C: bool = sys.argv[2].lower() == "t"

    generate_data(dimensions, C)


if __name__ == "__main__":
    main()
