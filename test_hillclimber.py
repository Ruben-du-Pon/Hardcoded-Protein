from code.algorithms.hillclimber import HillclimberFold
from code.classes.protein import Protein
from code.visualization import visualization_2D  # , visualization_3D

if __name__ == "__main__":

    test_protein = Protein("HHPHPHPHPHHHPPHHHHHHHHHH")
    # Create HillclimberFold object
    hillclimber = HillclimberFold(test_protein, 2, 100, True)

    # Run the hillclimber algorithm
    protein = hillclimber.run()

    protein.create_csv("data/output/csv/hillclimber_result.csv")
    visualization_2D.plot_2d(protein, ("red", "blue", "green"),
                             "data/output/plot/hillclimber_result.png", "png")
