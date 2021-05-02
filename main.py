
#----------------------
# Shay Fletcher     318727641
# Nika Tatikishvili 321735433
# https://github.com/shayfletcherz/NumericalAnalysisProject3.git
#----------------------

#יצירת מטריצת אפסים
def zeroMatrix(row, col):
    c = []
    for i in range(row):
        c += [[0] * col]
    return c

#כפל מטריצות
def multiplyMatrix(matrixA, matrixB):
    if len(matrixA[0]) is len(matrixB):
        c = zeroMatrix(len(matrixA), len(matrixB[0]))
        for row in range(len(matrixA)):
            for col in range(len(matrixB[0])):
                for x in range(len(matrixA)):
                    c[row][col] += (matrixA[row][x] * matrixB[x][col])
        return c
    return None

#פונקציה לכפל מטריצות בסקלר
def multScalar(matrix, scalarNum):
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            matrix[i][j] *= scalarNum
    return matrix

#חיבור מטריצות
def addMatrix(matrixA, matrixB):
    plusM = zeroMatrix(len(matrixA), len(matrixA[0]))
    for i in range(len(matrixA)):
        for j in range(len(matrixA[0])):
            plusM[i][j] = matrixA[i][j] + matrixB[i][j]
    return plusM

#הורדת שורה וטור
def remove(matrix, row, col):
    if row >= len(matrix) and col >= len(matrix):
        return matrix
    c = zeroMatrix(len(matrix) - 1, len(matrix) - 1)
    x = 0
    y = 0
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            if i is not row and j is not col:
                c[x][y] = matrix[i][j]
                if y is len(c[0]) - 1:
                    x += 1
                    y = 0
                else:
                    y += 1
    return c

#מציאת דטרמיננטה
def det(matrix):
    if len(matrix) and len(matrix[0]) is 2:
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
    sum1 = 0
    for j in range(len(matrix[0])):
        if j % 2 == 0:
            multBy = 1
        else:
            multBy = -1
        sum1 += multBy * matrix[0][j] * det(remove(matrix, 0, j))
    return sum1


#יצירת מטריצה אלמנטרית
def elementalMatrics(matrix, row, col):
    c = zeroMatrix(len(matrix), len(matrix[0]))
    c = identityMatrix(c)
    c[row][col] = -1 * (matrix[row][col] / matrix[col][col])
    return c

#יצירת מטריצת יחידה
def identityMatrix(i):
    for x in range(len(i[0])):
        i[x][x] = 1
    return i

#פונקציה להחלפה בין שורות שתקח את השורה ותחליף בשורה שנרצה
def swapRow(a, firstRow, swapRow):
    c = zeroMatrix(len(a), len(a[0]))
    c = identityMatrix(c)
    c[firstRow] = a[swapRow]
    c[swapRow] = a[firstRow]
    return c

#פונקציה שתמצא את המטריצה העליונה ואת מטריצה L ההופכית על מנת שנוכל לשנות אותה בהיפוך המטריצות
def findUpper(matrix):
    invL = identityMatrix(zeroMatrix(len(matrix), len(matrix[0])))
    for row in range(len(matrix)):
        j = row + 1
        while j < len(matrix):
            if matrix[row][row] is 0:  #פיבוט
                k = j + 1
                while k < len(matrix):
                    if matrix[k][row] is not 0:
                        b = swapRow(matrix, row, k)
                        matrix = multiplyMatrix(b, matrix)
                        invL = multiplyMatrix(b, invL)
                        break
                    k += 1
            b = elementalMatrics(matrix, j, row)
            matrix = multiplyMatrix(b, matrix)
            invL = multiplyMatrix(b, invL)
            j += 1
    return matrix, invL

#פונקציה לבדיקה האם יש למטריצה אלכסון דומיננטי
def dominantDiagonal(matrix):
    for i in range(len(matrix)):
        sum1 = 0
        for j in range(len(matrix[0])):
            if i is not j:
                sum1 += abs(matrix[i][j])
        if sum1 > matrix[i][i]:
            return False
    return True

#פונקציה שתחזיר מטריצה עם אלכסון של 1 ותעדכן את המטריצה
def diagonalOneMatrix(matrix, matInverse):
    b = zeroMatrix(len(matrix), len(matrix[0]))
    b = identityMatrix(b)
    for i in range(len(matrix[0])):
        b = identityMatrix(b)
        b[i][i] = 1 / matrix[i][i]
        matrix = multiplyMatrix(b, matrix)
        matInverse = multiplyMatrix(b, matInverse)
    return matrix, matInverse

#מחזירה מטריצה הופכית
def inverse(matrix):
    if det(matrix) is 0: #בדיקה האם ניתן לבצע הפיכה
        return
    matrix, matInverse = findUpper(matrix)
    matrix, matInverse = diagonalOneMatrix(matrix, matInverse)
    size = len(matrix[0]) - 1
    while size > 0: #הפיכה
        for i in range(size):
            b = elementalMatrics(matrix, i, size)
            matrix = multiplyMatrix(b, matrix)
            matInverse = multiplyMatrix(b, matInverse)
        size -= 1
    return matInverse

#פונקציה שבודקת את מקסימום סכומי השורה ומציגה את השורה המקסימלית
def infNorm(matrix):
    norm = 0
    for i in range(len(matrix[0])):
        sumRow = 0
        for j in range(len(matrix)):
            sumRow += abs(matrix[i][j])
        norm = max(sumRow, norm)
    return norm

#פונקציה לבדיקת התכנסות בשביל תנאי לחישוב - אם הנורם קטן מ1
def checkConvergence(matrix):
    if infNorm(matrix) < 1:
        return True
    return False

#פונקציה שתפרק את המטריצה למטריצה עליונה, תחתונה ואלכסונית
def calcLDU(a):

    L = zeroMatrix(len(a), len(a[0]))
    D = zeroMatrix(len(a), len(a[0]))
    U = zeroMatrix(len(a), len(a[0]))
    for i in range(len(a)):
        for j in range(len(a[0])):
            if i is j:
                D[i][j] = a[i][j] #יצירת מטריצה אלכסונית
            elif i < j:
                U[i][j] = a[i][j] #יצירת מטריצה עליונה
            else:
                L[i][j] = a[i][j] #יצירת מטריצה תחתונה
    return L, D, U

#ביצוע מטריצות GH של יעקובי
def GHjacobi(matrix):
    L, D, U = calcLDU(matrix)
    invD = inverse(D)
    G = multScalar(multiplyMatrix(invD, addMatrix(L, U)), -1)
    return G, invD

#פונקציה להצגת איטרציות בשיטת יעקובי
def calcJaacobi(matrix, vec):
    G, H = GHjacobi(matrix)
    print("\nJaacobi Method:\n ")
    if checkConvergence(G) is False:
        print("\nThe system can not calculate.")
        return
    guess(G, H, vec)

#ביצוע מטריצות GH של גאוס
def GHgauss(matrix):
    L, D, U = calcLDU(matrix)
    invLminusD = inverse(addMatrix(L, D))
    invLplusD = inverse(addMatrix(L, D))
    G = multScalar(multiplyMatrix(invLminusD, U), -1)
    return G, invLplusD

#חישוב גאוס והצגה בעזרת GH ו-וקטור תוצאה והצגת האיטרציות
def guess(G, H, vec):
    e = 0.00001 #תנאי עצירה
    iteration = 0
    oldNum = zeroMatrix(len(vec), len(vec[0]))
    newNum = zeroMatrix(len(vec), len(vec[0]))
    flag = True
    print("\n")
    while abs(newNum[0][0] - oldNum[0][0]) > e or flag is True: #ריצה על האיטרציות ובדיקה
        flag = False
        print(oldNum)
        newNum = oldNum
        oldNum = addMatrix(multiplyMatrix(G, oldNum), multiplyMatrix(H, vec))
        iteration += 1
    print("\nTotal number of iterations: " + str(iteration))

#חישוב בעזרת גאוס והצגת איטרציות
def calcGaussSeidel(matrix, vec):
    G, H = GHgauss(matrix) #שימוש בGH של גאוס
    print("\nGauss-Seidel Method\n ")
    if checkConvergence(G) is False:
        print("\nThe system can not calculate.")
        return
    guess(G, H, vec)


#def main():
