import random
import csv


def generate_random_sequence_H_P(length):
    return "".join(random.choice(["H", "P"]) for _ in range(length))


def generate_random_sequence_H_P_C(length):
    return "".join(random.choice(["H", "P", "C"]) for _ in range(length))


start_range = 1
end_range = 100
num_numbers = 50

additional_sequences_H_P = [
    generate_random_sequence_H_P(length)
    for length in random.sample(range(start_range, end_range + 1), num_numbers)
]

additional_sequences_H_P_C = [
    generate_random_sequence_H_P_C(length)
    for length in random.sample(range(start_range, end_range + 1), num_numbers)
]

all_sequences_H_P = [
    "HHPHHHPHPHHHPH",
    "HPHPPHHPHPPHPHHPPHPH",
    "PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP",
    "HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH",
] + additional_sequences_H_P

all_sequences_H_P_C = [
    "PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP",
    "CPPCHPPCHPPCPPHHHHHHCCPCHPPCPCHPPHPC",
    "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH",
    "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH",
] + additional_sequences_H_P_C


# Creating the CSV files
with open("sequences_H_P.csv", mode="w", newline="") as file:
    writer = csv.writer(file)

    for sequence in all_sequences_H_P:
        writer.writerow([sequence])

with open("sequences_H_P_C.csv", mode="w", newline="") as file:
    writer = csv.writer(file)

    for sequence in all_sequences_H_P_C:
        writer.writerow([sequence])
