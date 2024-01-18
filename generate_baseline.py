import csv
import sys
from .code.algorithms.random import RandomFold
from .code.classes.protein import Protein
from .code.classes.amino_acid import AminoAcid


def generate_baseline(dimensions: int, C: bool) -> None:

    if C:
        filename: str = "data/input/sequences_H_P_C.csv"
    else:
        filename: str = "data/input/sequences_H_P.csv"

    with open(filename, "r") as file:
        reader = csv.reader(file)
        line_number: int = 0

        for row in reader:
            if not row:
                break

            sequence: str = row[0]
            test_protein: Protein = Protein(sequence)
            test = RandomFold(test_protein, dimensions, True)
            scores = []

            if C:
                outputfile = f"data/output/baseline/baseline_{dimensions}D_C.csv"
            else:
                outputfile = f"data/output/baseline/baseline_{dimensions}D.csv"

            for iteration in range(1000):
                test.run()
                scores.append(test_protein.get_score())

                with open(outputfile, "a") as file:
                    writer = csv.writer(file)
                    writer.writerow(
                        [iteration, str(test_protein), test_protein.get_score()])
                    file.close()

            with open(outputfile, "a") as file:
                writer = csv.writer(file)
                writer.writerow(
                    ["", "", f"Average: {sum(scores) / len(scores)}"])
                file.close()

            line_number += 1


def main():
    if len(sys.argv) != 3:
        print("Usage: python generate_baseline.py <dimensions> <t/f>")
        sys.exit(1)

    dimensions: int = int(sys.argv[1])
    C: bool = sys.argv[2].lower() == "t"

    generate_baseline(dimensions, C)


if __name__ == "__main__":
    main()
