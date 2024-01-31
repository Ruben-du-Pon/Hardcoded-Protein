import csv
import sys
from ..algorithms.random import RandomFold
from ..classes.protein import Protein


def generate_baseline(dimensions: int, C: bool) -> None:
    """
    Generate a baseline for a protein folding problem.

    Parameters
    ----------
    dimensions : int
        The number of dimensions in which the protein folding will be simulated.
    C : bool
        A boolean indicating whether the protein sequence includes cysteine (C) or not.

    """
    # Set the input and output filenames
    if C:
        filename: str = "data/input/sequences_H_P_C.csv"
        outputfile: str = f"data/output/baseline/{dimensions}D_C.csv"
    else:
        filename: str = "data/input/sequences_H_P.csv"
        outputfile: str = f"data/output/baseline/{dimensions}D.csv"

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

            # Run the algorithm 10**5 times and write the results to the output file
            for iteration in range(10**5):
                test_protein: Protein = Protein(sequence)
                test = RandomFold(test_protein, dimensions, True)
                test_protein = test.run()
                scores.append(test_protein.get_score())

                with open(outputfile, "a") as file:
                    writer = csv.writer(file)
                    writer.writerow(
                        [iteration, str(test_protein), test_protein.get_score()])

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
    The main function that checks command line arguments and calls the generate_baseline function.
    """
    if len(sys.argv) != 3:
        print("Usage: python generate_baseline.py <dimensions> <t/f>")
        sys.exit(1)

    dimensions: int = int(sys.argv[1])
    C: bool = sys.argv[2].lower() == "t"

    generate_baseline(dimensions, C)


if __name__ == "__main__":
    main()
