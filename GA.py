import numpy as np
import random
import matplotlib.pyplot as plt
import pandas as pd
import time
from Fitness import fitness_fun


def enforce_five_ones(individual):
    num_ones = np.sum(individual)
    while num_ones > 5:
        # 随机选择一个1并将其设置为0
        one_indices = np.where(individual == 1)[0]
        individual[np.random.choice(one_indices)] = 0
        num_ones -= 1
    while num_ones < 5:
        # 随机选择一个0并将其设置为1
        zero_indices = np.where(individual == 0)[0]
        individual[np.random.choice(zero_indices)] = 1
        num_ones += 1
    return individual


def ini_pop(POPULATION_SIZE, GENE_LENGTH, ACTIVE_GENES):
    population = []
    for i in range(GENE_LENGTH):
        individual = np.zeros(GENE_LENGTH, dtype=int)
        individual[i] = 1
        population.append(individual)

    return population


# 初始化种群
def initialize_population(population_size, gene_length, active_genes):
    population = []
    for _ in range(population_size):
        individual = np.zeros(gene_length, dtype=int)
        active_indices = random.sample(range(gene_length), active_genes)
        individual[active_indices] = 1
        population.append(individual)
    return population


# 选择操作（轮盘赌方法）
def roulette_wheel_selection(population, fitness_scores):
    total_fitness = sum(fitness_scores)
    selection_probs = [f / total_fitness for f in fitness_scores]
    selected_indices = np.random.choice(range(len(population)), size=len(population), p=selection_probs)
    return [population[i] for i in selected_indices]


# 交叉操作
def crossover(parent1, parent2, crossover_rate, gene_length):
    if random.random() < crossover_rate:
        crossover_point = random.randint(1, gene_length - 1)
        child1 = np.concatenate([parent1[:crossover_point], parent2[crossover_point:]])
        child2 = np.concatenate([parent2[:crossover_point], parent1[crossover_point:]])
        # 确保每个后代个体恰好有5个1
        child1 = enforce_five_ones(child1)
        child2 = enforce_five_ones(child2)
        return child1, child2
    return parent1, parent2


# 变异操作
def mutation(individual, mutation_rate, active_genes):
    for i in range(len(individual)):
        if random.random() < mutation_rate:
            individual[i] = 1 if individual[i] == 0 else 0
    # 确保个体恰好有5个1
    individual = enforce_five_ones(individual)
    return individual

























