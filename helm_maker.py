import csv
import random
import argparse

def generate_helm_peptide(length):
    aa_list = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", 
               "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
    sequence = ".".join(random.choice(aa_list) for _ in range(length))
    return f"PEPTIDE1{{{sequence}}}$$$$V2.0"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--length", type=int, default=10)
    parser.add_argument("-n", "--number", type=int, default=5)
    parser.add_argument("-o", "--output", type=str, default="peptides.csv")
    args = parser.parse_args()
    with open(args.output, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["ID", "HELM_Sequence"])
        
        for i in range(1, args.number + 1):
            # Using a simple ID format to avoid confusion
            seq_id = f"ID_{i}"
            writer.writerow([seq_id, generate_helm_peptide(args.length)])

    print(f"File created: {args.output}")

if __name__ == "__main__":
    main()