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

    def binaryCountDigits(self, n):
        num = (2 **n)
        binary = [[] for _ in range(0 , num)]

        for i in range(0, num):
            temp = i

            if temp == 0:
                binary[i] = [0 for _ in range(0, n)]
            else:
                while temp >=1:
                    binary[i].insert(0, temp %2)
                    temp = temp //2
            #Check if the binary is n bits long
            while len(binary[i]) <n:
                binary[i].insert(0, 0)

        return binary
    
    def calcDistance(self, gMatrix):
        binaryMatrix = self.binaryCountDigits(self.data_size)

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
    
    def sydroms(self):
        data = [1, 1, 1, 1]

        gMatrix = self.G_matrix

        hMatrix = self.H_matrix
        hMatrix = hMatrix.transpose()

        codeword = np.matmul(data, gMatrix)
        for i in range(len(codeword)):
            if codeword[i] %2 == 0:
                codeword[i] = 0
            else:
                codeword[i] = 1

        sydromns = [[0] for _ in range(7)]

        for i in range(len(codeword)):
            if codeword[i] == 0:
                past = 0
                codeword[i] = 1
            else:
                past = 1
                codeword[i] = 0

            sydromn = np.matmul(codeword, hMatrix)
            for k in range(len(sydromn)):
                if sydromn[k] %2 == 0:
                    sydromn[k] = 0
                else:
                    sydromn[k] = 1
            sydromns[i] = sydromn
            codeword[i] = past

        c = 0
        count = 0

        for i in range(len(sydromns)):
            for j in range(i +1, len(sydromns)):
                if (sydromns[i] != sydromns[j]).any():
                    c += 1

            if c == len(sydromns) -(i +1):
                count += 1
            c = 0

        return count

    def fitness(self):
        self.fitness_value = 0
        self.matrices()
        #does H matrix have a standard amount of ones?
        ones = 0
        for i in range(0, len(self.H_matrix)):
            for j in range(0, len(self.H_matrix[i])):
                if self.H_matrix [i] [j] == 1:
                    ones += 1 #max 21
        if ones >= 12: #1 à 63, fazer contador binário
            self.fitness_value += 1
            #is G matrix linear?
            isGLI = self.isLinearIndepent(self.G_matrix) #Alterar dimensões
            if isGLI:
                self.fitness_value += 1
                #is H*G^T matrix all zeroes?
                transposeG = self.G_matrix.transpose()
                mult = np.matmul(self.H_matrix, transposeG)

                zeros = 0
                for i in range(0, len(mult)):
                    for j in range(0, len(mult[i])):
                        #garante apenas binarios
                        if mult [i] [j] %2 == 0:
                            mult [i] [j] = 0
                        if mult [i] [j] %2 == 1:
                            mult [i] [j] = 1
                        #Contar zeros
                        if mult [i] [j] == 0:
                            zeros += 1 #max 12
                if zeros == 12: #dados x redundancia
                    self.fitness_value += 12
                    #Minimun Distance
                    self.minDistanceValue = self.calcDistance(self.G_matrix)
                    self.fitness_value += self.minDistanceValue

                    if self.minDistanceValue > 2:
                        #There are 7 different sydroms?
                        if self.sydroms() == 7: #2^(n) -1
                            self.fitness_value += 1
                else:
                    self.fitness_value += zeros

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

        while f < 18:
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
            if population[i].fitness_value >= 18: 
                j += 1
                print("Indivíduo: ",population[i].genes)
                print("Matriz H: ", population[i].H_matrix)
                print("Matriz G: ",population[i].G_matrix)
                print("Distancia Minima: ",population[i].minDistanceValue)
                print("Fitness: ",population[i].fitness_value, "\n")
        return population

codeword_size = 7 #20
data_size = 4 #8

population_size = 1000
gens = 1000
mutation = 10 #Em porcentagem, chance de um idivíduo sofrer mutação
mutationH = 2 #Quantiade de bits
mutationG = 3
crossover = 50 

#Points table:
#1 - H matrix have a standard amount of ones or more (>= 12)
#2 - G matrix is linear
#2 -> 13 - H*G^T matrix != 0 
#14 - H*G^T matrix == 0 
#14 - Min distance = 0
#15 - Min distance = 1
#16 - Min distance = 2
#17 - Min distance = 3
#18 - Min distance = 3 +sydromns

main = Main(data_size, codeword_size, population_size, gens, mutation, crossover, mutationH, mutationG)
main.start()

#Matriz G 16x40
#Matriz H 24x40
