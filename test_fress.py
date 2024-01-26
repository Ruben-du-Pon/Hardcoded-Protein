from codefiles.algorithms.fress import FressFold
from codefiles.classes.protein import Protein
from codefiles.visualization import visualization_2D
from codefiles.visualization import visualization_3D
import time

if __name__ == "__main__":
    dimensions = 2
    test_protein = Protein(
        "HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH")
    # Create FressFold object
    fress = FressFold(test_protein, dimensions, 100, True)

    # Run the fress algorithm
    start_time = time.time()
    test_protein = fress.run() 
    end_time = time.time()  # Record the end time
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time} seconds")


    # Set the output filenames
    fold_algorithm = "FRESS"
    line_number = 0
    if dimensions == 2:
        filename = f"data/output/csv/{fold_algorithm}_{line_number}_2D.csv"
        plotname = f"data/output/plot/{fold_algorithm}_{line_number}_2D.png"
    else:
        filename = f"data/output/csv/{fold_algorithm}_{line_number}_3D.csv"
        plotname = f"data/output/plot/{fold_algorithm}_{line_number}_3D.png"

    test_protein.create_csv(filename)
    if dimensions == 2:
        visualization_2D.plot_2d(
            test_protein, ("red", "blue", "green"), plotname, "png")
    else:
        visualization_3D.plot_3d(
            test_protein, ("red", "blue", "green"), plotname, "png")
