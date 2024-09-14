import random
import sympy
import numpy as np

class geneticOperations:
    def mutation(self, population, mutationRate, mutatedBits):
        for i in range(len(population)):
            if random.randint(1, 100) <= mutationRate:
                k = 0
                while k < mutatedBits:
                    k += 1
                    j = random.randint(0, len(population[i].genes) -1)
                    if population[i].genes[j]:
                        population[i].genes[j] = 0
                    else:
                        population[i].genes[j] = 1
        return population

    def crossover(self, population, set):
        i = -1
        while len(population) < popSize//2:
            i += 1
            parent1 = population[i]
            #print("Pai 1: ")
            #print(parent1.genes)
            i += 1

            parent2 = population[i]
            #print("Pai 2: ")
            #print(parent2.genes)

            for _ in range(2):
                if set:
                    child_genome = gMatrix(dataSize, codewordSize)
                if set == 0:
                    child_genome = hMatrix(dataSize, codewordSize)

                k = 0
                while k < (parent1.size) *(cutPoint/100):
                    child_genome.genes[k] = parent1.genes[k]
                    k += 1
                while k < (parent2.size):
                    child_genome.genes[k] = parent2.genes[k]
                    k += 1
                #print("Filho: ")
                #print(child_genome.genes)
                population.append(child_genome)

                parent1, parent2 = parent2, parent1

        return population

def binaryCounter(num):
        binary = []
        while num >=1:
            binary.insert(0, num %2)
            num = num //2

        return binary

class gMatrix:
    def __init__(self, dataSize, codewordSize, fitnessValue = 0):
        self.dataSize = dataSize
        self.codewordSize = codewordSize
        self.size = self.dataSize *self.codewordSize
        self.fitnessValue = fitnessValue

        self.genes = [random.randint(0, 1) for _ in range(self.size)]

    def matrix(self):
        G_matrix = []

        for i in range(0, self.size):
            G_matrix.insert(i, self.genes[i])

        G_matrix = np.array(G_matrix)
        G_matrix = G_matrix.reshape((self.dataSize, self.codewordSize))

        return G_matrix

    def isLinearIndepent(self):
        matrix = self.matrix()
        qtd_linhas_independente = 0
        total_linhas = len(matrix)

        _, inds = sympy.Matrix(matrix).T.rref()
        for _ in inds:
            qtd_linhas_independente += 1

        if qtd_linhas_independente == total_linhas:
            return True
        else:
            return False

    def minimunDistance(self):
        binaryMatrix = []
        #Calculate all the poxible binary numbers
        for i in range(0, 2** self.dataSize):
            binaryMatrix.append(binaryCounter(i))
            #Make them all the same width
            while len(binaryMatrix[i]) < self.dataSize:
               binaryMatrix[i].insert(0, 0)
        
        gMatrix = self.matrix()
        gMatrix = np.array(gMatrix)
        
        codewords = [0 for _ in range(len(binaryMatrix))]
        for i in range(len(binaryMatrix)):
            codewords[i] = np.matmul(binaryMatrix[i], gMatrix)
        
        codewords = np.array(codewords)

        for i in range(len(codewords)):
            for j in range(7):
                codewords[i][j] = codewords[i][j] % 2
        
        distance = 0
        minDistance = 8
        for i in range(0, len(codewords)):
            if i == 0:
                zeros = 0
                for k in range(0, len(codewords[i])):
                    if codewords[i][k] == 0:
                        zeros += 1
                if zeros == 7:
                    i += 1

            for j in range(i + 1, len(codewords)):
                if j == len(codewords):
                    break
                for k in range(0, len(codewords[i])):
                    if codewords[i][k] != codewords[j][k]:
                        distance += 1
                if distance < minDistance:
                    minDistance = distance
                distance = 0
        return minDistance
    
    def fitness(self):
        self.fitnessValue = 0

        matrix = self.matrix()
        matrix = np.array(matrix)

        row = len(matrix)
        col = len(matrix[0])

        flag = True
        #Check if there is any column that is all 0s in G
        for j in range(col):
            if all(matrix[i][j] == 0 for i in range(row)):
                flag = False #False means that there is a column that is all 0s

        if flag:
            self.fitnessValue += 1

            if self.isLinearIndepent():
                self.fitnessValue += 1
                self.fitnessValue += self.minimunDistance()       

class hMatrix:
    def __init__(self, dataSize, codewordSize, fitnessValue = 0):
        self.dataSize = dataSize
        self.codewordSize = codewordSize
        self.redundancySize = codewordSize -dataSize
        self.size = self.redundancySize *self.codewordSize
        self.fitnessValue = fitnessValue

        self.genes = [random.randint(0, 1) for _ in range(self.size)]

    def matrix(self):
        H_matrix = []

        for i in range(0, self.size):
            H_matrix.insert(i, self.genes[i])

        H_matrix = np.array(H_matrix)
        H_matrix = H_matrix.reshape((self.redundancySize, self.codewordSize))

        return H_matrix

    def sydroms(self):
        data = binaryCounter((2 **self.dataSize) -1)

        gMatrix = gMat.matrix()

        hMatrix = self.matrix()
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
        self.fitnessValue = 0

        matrix = self.matrix()
        matrix = np.array(matrix)

        row = len(matrix)
        col = len(matrix[0])

        flag = True
        #Check if there is any column that is all 0s in H
        for j in range(col):
            if all(matrix[i][j] == 0 for i in range(row)):
                flag = False 

        if flag:
            self.fitnessValue += 1
            #Check if there is any row that is all 0s in H
            for row in matrix:
                if all(elem == 0 for elem in row):
                    flag = False 

            if flag:
                self.fitnessValue += 1
                #Check if any columns are equal
                for i in range(col):
                    for j in range(i+1, col):
                        if all(matrix[row][i] == matrix[row][j] for row in range(len(matrix))):
                            flag = False
                if flag:
                    self.fitnessValue += 1

                    identity = []
                    #Calculate the identity matrix
                    for i in range(len(matrix)):
                        binary = binaryCounter(2 **i)
                        #Makes every row the same width
                        while len(binary) < len(matrix):
                            binary.insert(0, 0)
                        identity.append(binary)

                    identity = np.array(identity)
                    matrix = matrix.transpose()
                    #Check if the identity matrix in in H, c represents the amount of columns that H has of identity
                    c = 0
                    for i in range(len(identity)):
                        for j in range(len(matrix)):
                            if (identity[i] == matrix[j]).all():
                                c += 1
                                break

                    matrix = matrix.transpose()
                    self.fitnessValue += c
                    if c == len(identity):
                        g = gMat.matrix()
                        g = g.transpose()
                        #G*H^T
                        mult = np.matmul(matrix, g)

                        zeros = 0
                        for i in range(0, len(mult)):
                            for j in range(0, len(mult[i])):
                                #Guarantees only binary numbers
                                if mult[i][j] %2 == 0:
                                    mult[i][j] = 0
                                if mult[i][j] %2 == 1:
                                    mult[i][j] = 1
                                #Count zeroes
                                if mult[i][j] == 0:
                                    zeros += 1 

                        self.fitnessValue += zeros
                        if zeros == self.dataSize *self.redundancySize:
                            self.fitnessValue += self.sydroms()

def startG():
    population = []
    for i in range(popSize):
        population.append(gMatrix(dataSize, codewordSize))

    generation = 1
    fitness = 0

    while fitness < 5:
    #while generation <= gens:
        for i in range(0, popSize):
            population[i].fitness()

        #Metade dos melhores se reproduzem
        #population = sorted(population, key=lambda gMatrix: gMatrix.fitnessValue, reverse=True)[:popSize //2]
        #Top 25% se reproduzem, o resto é aleatório
        population = sorted(population, key=lambda gMatrix: gMatrix.fitnessValue, reverse=True)[:popSize //4]

        population = genOp.crossover(population, 1)
        population = genOp.mutation(population, mutationRate, mutatedBits)

        while len(population) < popSize:
            population.append(gMatrix(dataSize, codewordSize))

        print("Geração ", generation)
        generation += 1

        for i in range(len(population)):
            if population[i].fitnessValue > fitness:
                fitness = population[i].fitnessValue               
        print("Maior fitness: ", fitness, "\n")

    generation -= 1
    population = sorted(population, key=lambda gMatrix: gMatrix.fitnessValue, reverse=True)
    print(population[0].genes)
    return population[0], generation

def startH():
    population = []
    for i in range(popSize):
        population.append(hMatrix(dataSize, codewordSize))

    generation = 1
    fitness = 0

    while fitness < 25:
    #while generation <= gens:
        for i in range(0, popSize):
            population[i].fitness()

        #Metade se reproduz
        #population = sorted(population, key=lambda hMatrix: hMatrix.fitnessValue, reverse=True)[:popSize //2]
        #Top 25% se reproduz
        population = sorted(population, key=lambda hMatrix: hMatrix.fitnessValue, reverse=True)[:popSize //4]

        population = genOp.crossover(population, 0)
        population = genOp.mutation(population, mutationRate, mutatedBits)

        while len(population) < popSize:
            population.append(hMatrix(dataSize, codewordSize))

        print("Geração ", generation)
        generation += 1

        for i in range(len(population)):
            if population[i].fitnessValue > fitness:
                fitness = population[i].fitnessValue               
        print("Maior fitness: ", fitness, "\n")

    generation -= 1
    population = sorted(population, key=lambda hMatrix: hMatrix.fitnessValue, reverse=True)
    print(population[0].genes)
    return population[0], generation

dataSize = 4
codewordSize = 7
popSize = 500

cutPoint = 80
mutationRate = 5
mutatedBits = 1

genOp = geneticOperations()
gens = 10

gensG = 0
gensH = 0

for _ in range(10):
    gMat, a = startG()
    gensG += a
    #gMat = gMatrix(dataSize, codewordSize)
    #gMat.genes = [1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1]

    hMat, b = startH()
    gensH += b

print("Matriz G: ", gensG/10)
print("Matriz H: ", gensH/10)

#[1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1]

#G points
#1 - No columns are all zero
#2 - Is linear independent
#5 - Minimun distance = 3

#H points
#1 - No columns are all zero
#2 - No rows are all zero
#3 - All columns are different of each other
#6 - Identity matrix is in H
#18 - H*G^T is all 0
#25 - The code have 7 distinct sydromns

#1000, 50, 10%
#Matriz G:  1.2
#Matriz H:  29.1

#1000, 50, 5%
#Matriz G:  1.1
#Matriz H:  22.9

#1000, 50, 1%
#Matriz G: 1.1
#Matriz H: 19.5

#500, 50, 10%
#Matriz G:  1.6
#Matriz H:  46.2

#500, 50, 5%
#Matriz G:  1.5
#Matriz H:  42.4

#500, 50, 1%
#Matriz G:  1.8
#Matriz H:  41.4

#1000, 80, 10%
#Matriz G:  1.1
#Matriz H:  19.5

#1000, 80, 5%
#Matriz G: 1.1
#Matriz H: 18.5

#1000, 80, 1%
#Matriz G: 1.1
#Matriz H: 18.9

#500, 80, 10%
#Matriz G: 1.4
#Matriz H: 42.2

#500, 80, 5%
#Matriz G: 1.3
#Matriz H: 36.4

#500, 80, 1%
#Matriz G: 1.6
#Matriz H: 44.5