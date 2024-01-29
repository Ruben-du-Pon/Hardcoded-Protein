from codefiles.classes.protein import Protein
from codefiles.algorithms.hillclimber import HillclimberFold
from codefiles.visualization import visualization_2D

if __name__ == "__main__":

    test_protein = Protein("PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP")
    # Create HillclimberFold object
    hillclimber = HillclimberFold(test_protein, 2, 500, verbose=True)

    # Run the hillclimber algorithm
    protein = hillclimber.run()
    protein.reset_grid()

    protein.create_csv("data/output/csv/hillclimber_data/result.csv")
    visualization_2D.plot_2d(protein, ("red", "blue", "green"),
                             "data/output/hillclimber_data/plots/result.png", "png")
