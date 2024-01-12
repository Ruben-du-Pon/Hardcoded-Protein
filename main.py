import csv
from code.classes import protein

def main():
    with open("data/input/sequences_H_P.csv", "r") as file:
        reader = csv.reader(file)
        line_number = 0
        
        for row in reader:
            if not row:
                break
            
            sequence= row[0]
            test_protein=protein.Protein(sequence)
            test_protein.create_csv(line_number)
            line_number += 1
            
if __name__ == "__main__":
    main()