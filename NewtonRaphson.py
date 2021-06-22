# Cálculo da raiz quadrada por newton raphson
# definir a funcao: no caso será f(x) = raiz quadrada de x
# definir um x0 (nossa escolha) e numero de iteracoes

# queremos saber se é mais rapido e mais preciso o calculo de 1/raiz de A
# se calcular primeiro raiz de A pelo newton raphson e depois inverter
# ou calcular a funcao 1/raiz de A direto pelo newton raphson


import numpy as np
import sys
import math
import struct

getBin = lambda x: x > 0 and str(bin(x))[2:] or "-" + str(bin(x))[3:]


def floatToBinary64(value):
    val = struct.unpack('Q', struct.pack('d', value))[0]
    return getBin(val)


def padraoIEEE(x):
    binstr = floatToBinary64(x)

    numerodeZeros = 64 - len(binstr)

    aux = ''

    for i in range(numerodeZeros):
        aux += '0'

    for i in binstr:
        aux += i

    binstr = aux

    S = int(binstr[0])
    E = binstr[1:12]
    M = binstr[12:]

    eResult = 0
    j = 0

    for i in range(10, -1, -1):
        eResult += int(E[j]) * pow(2, i)
        j += 1

    j = 0
    mResult = 0
    for i in range(51, -1, -1):
        mResult += int(M[j]) * pow(2, i)
        j += 1

    y = mResult / pow(2, 52)

    x0 = 1 + y/2

    ezinho = eResult - 1023

    return x0, ezinho, y
    
 
def newtonRaphson(A, x0):
    A = float(A)

    precisao = 5e-12
    erro = sys.maxsize
    itmax = 5
    
    k = 0
    xk = float(x0)

    resultado_esperado = math.sqrt(A)

    print("Iterações do NR-Convencional")

    while (erro >= precisao) and (k <= itmax):
        xk1 = 1/2 * (xk + A/xk)
        erro = abs(xk - xk1)
        
        print("\tk = ", str(k), "  xk1 = ", str(xk1), "  desvio = ", str(abs(resultado_esperado - xk1)), "  errou", erro)

        xk = xk1
        k = k + 1

    return xk


def newtonRaphsonInverso(A,x0):
    A = float(A)

    precisao = 5e-12
    erro = sys.maxsize
    itmax = 5
    
    k = 0
    xk = float(x0)

    resultado_esperado = 1/math.sqrt(A)

    print("Iterações do NR-Inverso")

    while (erro >= precisao) and (k <= itmax):
        xk1 = abs(xk * (1.5 - ((0.5 * A) * (xk * xk))))
        erro = abs(xk - xk1)
        
        print("\tk = ", str(k), "  xk1 = ", str(xk1), "  desvio = ", str(abs(resultado_esperado - xk1)), "  errou = ", erro)

        xk = xk1
        k = k + 1

    return xk



def main():
    radicando = input("Digite o radicando: ")

    x0, expoente, f = padraoIEEE(float(radicando))

    resultadoConvencional = 0
    resultadoInverso = 0
    raizDe2 = 1.41421356237309504880168872420969807856967187537694807317667973799

    if(expoente % 2 == 0):
        resultadoConvencional = (pow(2, expoente/2)) * newtonRaphson(1 + f, x0)
        resultadoInverso = newtonRaphsonInverso(1 + f, x0) / (pow(2, expoente/2))

    else:
        resultadoConvencional = (pow(2, ((expoente -1)/2))) * raizDe2 * newtonRaphson(1 + f, x0)
        resultadoInverso = newtonRaphsonInverso(1 + f, 1/x0) / ((pow(2, ((expoente -1)/2))) * raizDe2)

    print("Resultado NR-Convencional: ", 1/resultadoConvencional)
    print("Resultado NR-Inverso: ", resultadoInverso)

main()
