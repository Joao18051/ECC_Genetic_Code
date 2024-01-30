#Test to check the minimum distance of the codewords
import numpy as np

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

gMatrix = [[0, 0, 1, 0, 1, 1, 1],
 [1, 0, 0, 0, 1, 1, 0],
 [0, 1, 1, 1, 0, 0, 1],
 [0, 1, 0, 0, 0, 1, 1]]

gMatrix = np.array(gMatrix)

codewords = [0 for _ in range(len(binaryMatrix))]

for i in range(len(binaryMatrix)):
    codewords[i] = np.matmul(binaryMatrix[i], gMatrix)

codewords = np.array(codewords)

for i in range(len(codewords)):
    for j in range(7):
        codewords[i] [j] = codewords[i] [j] %2

print(codewords, "\n")

distance = 0
minDistance = 8
linha1 = 0 #Those two are just to say wich rows have the said min distance
linha2 = 0 #It returns the last two rows that it found with the said distance
for i in range(0, len(codewords)):
    if i == 0:
        zeros = 0
        for k in range(0, len(codewords[i])):
            if codewords [i] [k] == 0:
                zeros += 1
        if zeros == 7:
            i += 1

    for j in range(i +1, len(codewords)):
        if j == len(codewords):
            break;
        for k in range(0, len(codewords[i])):
            if codewords[i] [k] != codewords[j] [k]:
                distance += 1
        if distance < minDistance:
            minDistance = distance
            linha1 = i
            linha2 = j
        distance = 0

print("Distancia minima: ", minDistance)
print("Linha 1: ", linha1)
print("Linha 2: ", linha2)