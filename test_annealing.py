from codefiles.classes.protein import Protein
from codefiles.algorithms.annealing import AnnealingFold
from codefiles.visualization import visualization_2D

if __name__ == "__main__":

    test_protein = Protein("PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP")
    # Create HillclimberFold object
    annealing = AnnealingFold(test_protein, 2, 250, verbose=True)

    # Run the hillclimber algorithm
    protein = annealing.run()
    protein.reset_grid()

    protein.create_csv("data/output/annealing_data/result.csv")
    visualization_2D.plot_2d(protein, ("red", "blue", "green"),
                             "data/output/annealing_data/plot/result.png", "png")
