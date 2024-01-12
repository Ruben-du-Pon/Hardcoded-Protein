from code.classes import protein

sequence = input("Enter the sequence of Aminoacids for the protein: ")
test_protein = protein.Protein(sequence)
test_protein.create_csv