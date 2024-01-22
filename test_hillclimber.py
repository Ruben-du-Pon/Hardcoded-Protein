from code.algorithms.hillclimber import HillclimberFold
from code.classes.protein import Protein

if __name__ == "__main__":

    test_protein = Protein("HHPHPHPHPHHHPPHHHHHHHHHH")
    # Create HillclimberFold object
    hillclimber = HillclimberFold(test_protein, 2, 100, True)

    # Run the hillclimber algorithm
    protein = hillclimber.run()

    protein.create_csv("data/output/csv/test_hillclimber.csv")
