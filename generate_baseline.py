import csv
from code.algorithms.random import RandomFold
from code.classes.protein import Protein
from code.visualization import visualization_2D


highscore = (None, 0)
highscore_plotname = None
valid_proteins = []


def experiment():
    global highscore
    global highscore_plotname
    global valid_proteins

    for iteration in range(50000):
        with open("data/input/sequences_H_P.csv", "r") as file:
            reader = csv.reader(file)
            line_number = 0

            for row in reader:
                if not row:
                    break

                sequence = row[0]
                protein = Protein(sequence)
                test = RandomFold(protein, 2, False)
                test.run()

                if check_valid(protein):
                    valid_proteins.append(protein)

                filename = f"data/output/baseline/csv/random_{line_number}_{iteration}.csv"
                plotname = f"data/output/baseline/plot/random_{line_number}_{iteration}.png"
                protein.create_csv(filename)
                visualization_2D.plot_2d(
                    protein, ("red", "blue", "green"), plotname)

                score = protein.get_score()
                if score > highscore[1]:
                    highscore = (protein, score)
                    highscore_plotname = plotname

                line_number += 1


def check_valid(protein):
    positions = set()
    current = protein.get_list()

    if current.position not in positions:
        positions.add(current.position)
    else:
        return False

    return True


def write_stats():
    global highscore
    global highscore_plotname
    global valid_proteins

    statsfile = "data/output/baseline/stats.md"
    with open(statsfile, "w") as file:
        file.write(f"# Baseline\n")
        file.write(f"Number of tests done: {50000}\n")
        file.write(f"Number of valid proteins: {len(valid_proteins)}\n")
        file.write(
            f"Highest scoring protein: {str(highscore[0])}, {highscore[1]}\n")
        file.write(f"![Highscore]({highscore_plotname})\n")


if __name__ == "__main__":
    experiment()
    write_stats()
