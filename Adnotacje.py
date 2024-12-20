import warnings
from Bio import SeqIO, BiopythonParserWarning
import re

warnings.simplefilter("ignore", BiopythonParserWarning)

def fix_locus(input_file, output_file):
    pattern = r"LOCUS\s+(\S+)_length_(\d+)_cov_\S+\s+\d*\s*bp.*"
    with open(input_file, "r") as file, open(output_file, "w") as outfile:
        for line in file:
            if line.startswith("LOCUS"):
                match = re.match(pattern, line)
                if match:
                    sequence_name = match.group(1)
                    sequence_length = match.group(2)
                    fixed_line = f"LOCUS       {sequence_name:<30} {sequence_length:>7} bp    DNA     linear\n" #name od 13 kolumny
                    outfile.write(fixed_line)
                    print(f"Naprawione LOCUS: {fixed_line.strip()}")
                else:
                    print(f"Nie rozpoznano LOCUS: {line.strip()}")
            else:
                outfile.write(line)

def analyze_genes(file_path):
    count_genes = 0
    total_gene_length = 0
    polymerase_genes = 0

    for record in SeqIO.parse(file_path, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":  
                count_genes += 1
                gene_length = len(feature.location)
                total_gene_length += gene_length
                if "product" in feature.qualifiers:
                    product = feature.qualifiers["product"][0].lower()
                    if "polymerase" in product:
                        polymerase_genes += 1
    avg_length = total_gene_length/count_genes 
    return count_genes, avg_length, polymerase_genes

"""
def extract_polymerase_file(file_path, output_file):
    polymerase_genes = []
    for record in SeqIO.parse(file_path, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                if "product" in feature.qualifiers:
                    product = feature.qualifiers["product"][0].lower()
                    if "polymerase" in product:
                        polymerase_genes.append(feature)
    with open(output_file, "w") as file:
        for gene in polymerase_genes:
            file.write(str(gene) + "\n")
"""


def compare_annotations(file_bakta, file_prokka):
    fixed_prokka_file = "fixed_prokka.gbk"
    fix_locus(file_prokka, fixed_prokka_file)


    print("Analiza pliku BAKTA:")
    count_bakta, avg_length_bakta, polymerase_bakta = analyze_genes(file_bakta)
    #extract_polymerase_file(file_bakta, "bakta_polymerases.txt")
    print(f"Liczba genów: {count_bakta}")
    print(f"Średnia długość genów: {avg_length_bakta:.2f}")
    print(f"Liczba genów polimerazy: {polymerase_bakta}")
    print("-" * 50)

    # Analiza i zapis genów dla PROKKA
    print("Analiza pliku PROKKA:")
    count_prokka, avg_length_prokka, polymerase_prokka = analyze_genes(fixed_prokka_file)
    #extract_polymerase_file(fixed_prokka_file, "prokka_polymerases.txt")
    print(f"Liczba genów: {count_prokka}")
    print(f"Średnia długość genów: {avg_length_prokka:.2f}")
    print(f"Liczba genów polimerazy: {polymerase_prokka}")


file_bakta = "C:/Users/DELL 5470/Downloads/assembly.gbff"
file_prokka = "C:/Users/DELL 5470/Downloads/assembly.gbk"

compare_annotations(file_bakta, file_prokka)
