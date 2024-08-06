from Bio import Entrez, SeqIO, AlignIO
from Bio.Align import substitution_matrices
import numpy as np
import random
import matplotlib.pyplot as plt


# Data retrieval and handling #
def fetch_sequences_from_genbank(accession_numbers):
    Entrez.email = "fredjjones1700@gmail.com"
    sequences = []
    for acc in accession_numbers:
        handle = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text")
        seq_record = SeqIO.read(handle, "fasta")
        sequences.append(str(seq_record.seq))
        handle.close()
    return sequences

# Example accession numbers
accession_numbers = ["NM_001301716" ,"NM_001301717", "NM_001301718", "NM_001301719", "NM_001301720", "NM_001301721"]
sequences = fetch_sequences_from_genbank(accession_numbers)

##########################################################################################
# Genetic Algorithm

# Generate initial population
def generate_initial_population(sequences, population_size):
    max_length = max(len(seq) for seq in sequences)
    population = []
    for _ in range(population_size):
        alignment = [list(seq.ljust(max_length, '-')) for seq in sequences]
        # Add random gaps
        for seq in alignment:
            num_gaps = random.randint(0, len(seq) // 2)
            for _ in range(num_gaps):
                position = random.randint(0, len(seq))
                seq.insert(position, '-')
        population.append(alignment)
    return population

# Fitness function using BLOSUM62 substitution matrix
def calculate_fitness(alignment):
    substitution_matrix = substitution_matrices.load('BLOSUM62')
    score = 0
    alignment_length = len(alignment[0])
    
    for i in range(alignment_length):
        column = [seq[i] for seq in alignment if i < len(seq)]
        for j in range(len(column) - 1):
            for k in range(j + 1, len(column)):
                if column[j] != '-' and column[k] != '-':
                    score += substitution_matrix[column[j], column[k]]
                elif column[j] == '-' or column[k] == '-':
                    score -= 1  # Penalty for gaps, adjust this value as needed
    return score

# Selection
def roulette_wheel_selection(population, fitnesses):
    min_fitness = min(fitnesses)
    shifted_fitnesses = [f - min_fitness for f in fitnesses]
    total_fitness = sum(shifted_fitnesses)
    
    if total_fitness == 0:
        # Assign equal probability if total_fitness is zero
        probabilities = [1 / len(population)] * len(population)
    else:
        probabilities = [f / total_fitness for f in shifted_fitnesses]
        
    selected_indices = np.random.choice(len(population), size=len(population), p=probabilities)
    return [population[i] for i in selected_indices]

# Crossover
def crossover(parent1, parent2):
    point = random.randint(1, len(parent1[0]) - 2)
    child1 = [p1[:point] + p2[point:] for p1, p2 in zip(parent1, parent2)]
    child2 = [p2[:point] + p1[point:] for p1, p2 in zip(parent1, parent2)]
    return child1, child2

# Mutation
def mutate(alignment, mutation_rate):
    for seq in alignment:
        if random.random() < mutation_rate:
            position = random.randint(0, len(seq) - 1)
            seq[position] = random.choice("ACGT-")
    return alignment

# Replacement
def replace_population(old_population, new_population):
    return new_population

def elite_selection(population, fitnesses, elite_size):
    sorted_population = [x for _, x in sorted(zip(fitnesses, population), key=lambda pair: pair[0], reverse=True)]
    return sorted_population[:elite_size]

# Main GA loop
def run_genetic_algorithm(sequences, population_size=100, num_generations=1000, mutation_rate=0.1, elite_size=5):
    # Generate initial population
    population = generate_initial_population(sequences, population_size)
    best_fitness_over_time = []
    previous_best_fitness = None
    
    for generation in range(num_generations):
        fitnesses = [calculate_fitness(alignment) for alignment in population]
        current_best_fitness = max(fitnesses)
        best_fitness_over_time.append(current_best_fitness)
        
        
        previous_best_fitness = current_best_fitness
        elite = elite_selection(population, fitnesses, elite_size)
        selected_population = roulette_wheel_selection(population, fitnesses)
        new_population = elite
        
        for i in range(0, len(selected_population) - elite_size, 2):
            parent1, parent2 = selected_population[i], selected_population[i + 1]
            child1, child2 = crossover(parent1, parent2)
            new_population.extend([mutate(child1, mutation_rate), mutate(child2, mutation_rate)])
        
        population = replace_population(population, new_population)
    
    best_alignment = max(population, key=calculate_fitness)
    return best_alignment, best_fitness_over_time

# Function to measure column conservation
def calculate_column_conservation(alignment):
    alignment_length = len(alignment[0])
    conservation_scores = []
    
    for i in range(alignment_length):
        column = [seq[i] for seq in alignment if i < len(seq)]
        unique_elements = set(column)
        conservation_score = len(unique_elements) / len(column)
        conservation_scores.append(conservation_score)
    
    return conservation_scores

# Function to visualize the alignments
def plot_fitnesses(fitness_over_time):
    generations = list(range(len(fitness_over_time)))
    
    plt.figure(figsize=(14, 10))
    plt.plot(generations, fitness_over_time, marker='o', linestyle='-', color='b')
    plt.title('Fitness Over Generations')
    plt.xlabel('Generation')
    plt.ylabel('Fitness')
    plt.grid(True)
    plt.show()

# Evaluate the GA performance
def evaluate_ga_performance(best_alignment, fitnesses):
    # Calculate the total alignment score
    total_score = calculate_fitness(best_alignment)
    print("Total Alignment Score:", total_score)
    
    
    # Visualize the alignment
    plot_fitnesses(fitnesses)

# Read ClustalW alignment
def read_clustalw_alignment(output_file):
    alignment = AlignIO.read(output_file, "clustal")
    aligned_sequences = [str(record.seq) for record in alignment]
    return aligned_sequences

# Convert alignment to list format
def convert_alignment_to_list(alignment):
    return [list(seq) for seq in alignment]

# Path to your ClustalW output file
clustalw_output_file = "clustal.fasta"

# Read and convert ClustalW alignment
clustalw_alignment = read_clustalw_alignment(clustalw_output_file)
clustalw_alignment_list = convert_alignment_to_list(clustalw_alignment)

# Calculate the SP score for ClustalW alignment
clustalw_sp_score = calculate_fitness(clustalw_alignment_list)
print(f"SP Score for ClustalW Alignment: {clustalw_sp_score}")


# Fetch sequences from GenBank
sequences = fetch_sequences_from_genbank(accession_numbers)

# Run GA for MSA with fetched sequences and evaluate performance
best_alignment, fitnesses = run_genetic_algorithm(sequences, population_size=200, num_generations=2000, mutation_rate=0.1)
best_alignment_list = ["".join(seq) for seq in best_alignment]


# Evaluate GA performance
evaluate_ga_performance(best_alignment, fitnesses)