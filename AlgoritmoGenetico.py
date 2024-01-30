import random
import sympy
import numpy as np

class Genome:
    def __init__(self, data_size, codeword_size, mutationH, mutationG, fitness_value=0):
        self.fitness_value = fitness_value
        self.minDistanceValue = 0

        self.codeword_size = codeword_size
        self.data_size = data_size
        self.parity_size = codeword_size -data_size

        self.H_size = codeword_size *self.parity_size
        self.G_size = codeword_size *data_size
        self.H_matrix = []
        self.G_matrix = []
        self.mutationH = mutationH
        self.mutationG = mutationG

        self.genes = [random.randint(0, 1) for _ in range(self.H_size +self.G_size)]
        
    def matrices(self):
        H_matrix = []
        G_matrix = []
        #Referente a matriz H
        for i in range(0, self.H_size):
            H_matrix.insert(i, self.genes[i])

        self.H_matrix = np.array(H_matrix)
        self.H_matrix = self.H_matrix.reshape((self.parity_size, self.codeword_size))

        #Referente a matriz G
        for j in range(0, self.G_size):
            G_matrix.insert(j, self.genes[i +j])

        self.G_matrix = np.array(G_matrix)
        self.G_matrix = self.G_matrix.reshape((self.data_size, self.codeword_size))
    
    def isLinearIndepent(self, matrix):
        qtd_linhas_independente = 0
        total_linhas = len(matrix)

        _, inds = sympy.Matrix(matrix).T.rref()
        for _ in inds:
            qtd_linhas_independente += 1

        if qtd_linhas_independente == total_linhas:
            return True
        else:
            return False

    def calcDistance(self, gMatrix):
        binaryMatrix = [[0, 0, 0, 0],
                [0, 0, 0, 1],
                [0, 0, 1, 0],
                [0, 0, 1, 1],
                [0, 1, 0, 0],
                [0, 1, 0, 1],
                [0, 1, 1, 0],
                [0, 1, 1, 1],
                [1, 0, 0, 0],
                [1, 0, 0, 1],
                [1, 0, 1, 0],
                [1, 0, 1, 1],
                [1, 1, 0, 0],
                [1, 1, 0, 1],
                [1, 1, 1, 0],
                [1, 1, 1, 1]
                ]

        gMatrix = np.array(gMatrix)

        codewords = [0 for _ in range(len(binaryMatrix))]

        for i in range(len(binaryMatrix)):
            codewords[i] = np.matmul(binaryMatrix[i], gMatrix)

        codewords = np.array(codewords)

        for i in range(len(codewords)):
            for j in range(7):
                codewords[i] [j] = codewords[i] [j] %2

        distance = 0
        minDistance = 8

        for i in range(0, len(codewords)):
            for j in range(i +1, len(codewords)):
                if j == len(codewords):
                    break;
                for k in range(0, len(codewords[i])):
                    if codewords[i] [k] != codewords[j] [k]:
                        distance += 1
                if distance < minDistance:
                    minDistance = distance
                distance = 0
        return minDistance
    
    def fitness(self):
        self.fitness_value = 0
        self.matrices()
        transposeG = self.G_matrix.transpose()
        mult = np.matmul(self.H_matrix, transposeG)

        zeros = 0
        for i in range(0, len(mult)):
            for j in range(0, len(mult[i])):
                if mult [i] [j] == 0:
                    zeros += 1 #max 12
        if zeros == 12:
            self.fitness_value += 20
        else:
            self.fitness_value += zeros 

        if zeros == 12:
            self.minDistanceValue = self.calcDistance(self.G_matrix)
            self.fitness_value += 3**self.minDistanceValue

        #isHLI = self.isLinearIndepent(self.H_matrix)
        #if isHLI:
            #self.fitness_value += 5

        isGLI = self.isLinearIndepent(self.G_matrix)
        if isGLI:
            self.fitness_value += 5

    def mutate(self, mutation_rate):
        if random.randint(1, 100) < mutation_rate:
            k = 0
            while k < mutationH:
                k += 1
                j = random.randint(0, self.H_size -1)
                if self.genes[j]:
                    self.genes[j] = 0
                else:
                    self.genes[j] = 1
            
            k = 0
            while k < mutationG:
                k += 1
                j = random.randint(self.H_size, self.H_size +self.G_size -1)
                if self.genes[j]:
                    self.genes[j] = 0
                else:
                    self.genes[j] = 1
class Main:
    def __init__(self, data_size, codeword_size, population_size, gens, mutation, crossover, mutationH, mutationG):
        self.codeword_size = codeword_size
        self.data_size = data_size

        self.population_size = population_size
        self.gens = gens

        self.mutation = mutation
        self.crossover = crossover

        self.mutationH = mutationH
        self.mutationG = mutationG

    def population(self):
        return [Genome(self.data_size, self.codeword_size, self.mutationH, self.mutationG) for _ in range(self.population_size)]
    
    def genarations(self):
        population = self.population()
        k = 1
        f = 0

        while f < 52:
        #for _ in range(self.gens):
            for i in range(self.population_size):
                population[i].fitness()
            population = sorted(population, key=lambda genome: genome.fitness_value, reverse=True)[:self.population_size //2]
            
            population = self.new_gen(population)#self.new_gen(newPop)
            print("Geração ", k)
            k += 1

            for i in range(self.population_size):
                if population[i].fitness_value > f:
                    f = population[i].fitness_value
            if f >= 28:
                print("Maior fitness: ", f)
        return population

    def new_gen(self, population: Genome):
        i = 0
        while len(population) < self.population_size:
            parent1 = population[i]
            i += 1

            if i > len(population[i].genes) - 1:
                i = 0
            parent2 = population[i]
            i += 1

            child_genome = Genome(data_size, codeword_size, mutationG, mutationH)
            k = 0
            while k < (parent1.H_size -1) *(self.crossover/100):
                child_genome.genes[k +7] = parent1.genes[k +7]
                k += 1
            if k != parent1.H_size:
                while k < (parent2.H_size -1) *(self.crossover/100):
                    child_genome.genes[k +7] = parent2.genes[k +7]
                    k += 1
            m = 0
            while m < (parent1.G_size -1) *(self.crossover/100):
                child_genome.genes[m +k +7] = parent1.genes[m +k +7]    
                m += 1
            if m != parent1.G_size:
                while m < (parent2.G_size -1) *(self.crossover/100):
                    child_genome.genes[m +k +7] = parent2.genes[m +k +7]    
                m += 1
                
            child_genome.mutate(self.mutation)
            population.append(child_genome)
            #population = np.append(population, child_genome)

        return population
    
    def start(self):
        population = self.genarations()

        j = 0
        for i in range(len(population) -1):
            if population[i].fitness_value >= 52: 
                j += 1
                print("Indivíduo: ",population[i].genes)
                print("Matriz H: ", population[i].H_matrix)
                print("Matriz G: ",population[i].G_matrix)
                print("Distancia Minima: ",population[i].minDistanceValue)
                print("Fitness: ",population[i].fitness_value, "\n")

        print(f'Há {j} combinações de matrizes que satisfazem as condições')

codeword_size = 7
data_size = 4

population_size = 1000
gens = 200
mutation = 10 #Em porcentagem, chance de um idivíduo sofrer mutação
mutationH = 2 #Quantiade de bits
mutationG = 1
crossover = 50 

#Tabela de pontuações:
#Se H *G^T não for uma matriz 0, fitness += qtd. de zeros
#Se H *G^T for uma matriz 0, fitness += 20
#Se G for L.I, fitness += 5
#A pontuação da distancia é dado por: fitness += 3^distanciaMinima
#(Distancia minima só é calculada se H *G^T = 0)
#Assim, um indivíduo que tem H *G^T = 0 e G L.I, tem fitness:
#26 - Distancia 0
#28 - Distancia 1
#34 - Distancia 2
#52 - Distancia 3
main = Main(data_size, codeword_size, population_size, gens, mutation, crossover, mutationH, mutationG)
main.start()
