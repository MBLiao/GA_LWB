import numpy as np
import random
import matplotlib.pyplot as plt
import pandas as pd
import time


def initialize_pop_loose(population_size, gene_length):
    population = []
    for _ in range(population_size):
        individual = np.zeros(gene_length, dtype=int)
        for i in range(gene_length):
            if random.random() < 0.1:
                individual[i] = 1
        population.append(individual)
    return population


# 选择操作（轮盘赌方法）
def roulette_wheel_selection(population, fitness_scores):
    total_fitness = sum(fitness_scores)
    selection_probs = [f / total_fitness for f in fitness_scores]
    selected_indices = np.random.choice(range(len(population)), size=20, p=selection_probs)
    return [population[i] for i in selected_indices]


# 交叉操作
def crossover(parent1, parent2, crossover_rate, gene_length):
    fg = 0
    if random.random() < crossover_rate:
        crossover_point = random.randint(1, gene_length - 1)
        child1 = np.concatenate([parent1[:crossover_point], parent2[crossover_point:]])
        child2 = np.concatenate([parent2[:crossover_point], parent1[crossover_point:]])
        fg = 1
        return child1, child2, fg
    return parent1, parent2, fg


# 变异操作
def mutation(individual, mutation_rate, gene_length, per_gene_variation=False):
    mflag = 0
    if per_gene_variation:
        for i in range(len(individual)):
            if random.random() < mutation_rate:
                individual[i] = 1 if individual[i] == 0 else 0
                mflag = 1
    else:
        if random.random() < mutation_rate:
            mp = random.randint(0, gene_length - 1)
            individual[mp] = 1 if individual[mp] == 0 else 0
            mflag = 1

    return individual, mflag


def select_fixed_size_population(population, fitness_scores, fixed_size):
    # 确保种群和适应度分数列表的长度相同
    assert len(population) == len(fitness_scores), "Population and fitness scores must be of the same length."
    population_with_scores = list(zip(population, fitness_scores))

    sorted_population_with_scores = sorted(population_with_scores, key=lambda x: x[1], reverse=True)
    selected_population_with_scores = sorted_population_with_scores[:fixed_size]
    selected_population, selected_fitness_scores = zip(*selected_population_with_scores)

    return list(selected_population), list(selected_fitness_scores)


def ini_rand_operation(population_size, gene_length):
    population = []
    for _ in range(population_size):
        individual = np.zeros(gene_length, dtype=int)
        for i in range(gene_length):
            o = random.random()
            if o < 0.03:
                individual[i] = 1
            elif o < 0.06:
                individual[i] = 2
            elif o < 0.09:
                individual[i] = 3
        population.append(individual)
    return population


def mutation_rand(individual, mutation_rate, gene_length, per_gene_variation=False):
    mflag = 0
    if per_gene_variation:
        for i in range(len(individual)):
            if random.random() < mutation_rate:
                individual[i] = 1 if individual[i] == 0 else 0
                mflag = 1
    else:
        if random.random() < mutation_rate:
            for i in range(len(individual)):
                om = random.random()
                if om < 0.01:
                    individual[i] = 1
                elif om < 0.02:
                    individual[i] = 2
                elif om < 0.03:
                    individual[i] = 3
                elif om < 0.05:
                    individual[i] = 0
            mflag = 1

    return individual, mflag















