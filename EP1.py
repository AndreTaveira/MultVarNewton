import math

def main():
  
    sair = False
    
    while(not sair):
        
        menu = int(input("Selecione vizualizar testes (1) ou localização por gps (2): "))
        
        if menu == 1:
            
            menu_teste = int(input("Selecione teste 1 (1), 2 (2) ou 3 (3): "))


            if menu_teste == 1:

                print("Chute inicial: [0, 0] \neps = 1e-5 \nmaximo de iterações: 50\n\n")
                k1, x1 = newton_multivar(Ft1, Jft1, [0, 0], 1e-5, 50)
                print("Número de iterações: %d" %k1)
                print("x:\n")
                print_vector(x1)
                print("F(x):\n")
                print_vector(Ft1(x1))
                
            if menu_teste == 2:

                print("Chute inicial: [1, 1, 1, 1] \neps = 1e-5 \nmaximo de iterações: 50\n\n")
                k2, x2 = newton_multivar(Ft2, Jft2, [1, 1, 1, 1], 1e-5, 50)
                print("Número de iterações: %d" %k2)
                print("x:\n")
                print_vector(x2)
                print("F(x):\n")
                print_vector(Ft2(x2))

            if menu_teste == 3:

                print("Chute inicial: vetor nulo de 20, 40 e 80 \neps = 1e-5 \nmaximo de iterações: 50")
                k31, x31 = newton_multivar(Ft3, Jft3, [1 for i in range(20)], 1e-5, 50)
                print("\n\nn = 20")
                print("Número de iterações: %d" %k31)
                print("x:")
                print_vector(x31)
                print("F(x):")
                print_vector(Ft3(x31))

                k32, x32 = newton_multivar(Ft3, Jft3, [1 for i in range(40)], 1e-5, 50)
                print("\n\nn = 40")
                print("Número de iterações: %d" %k32)
                print("x:")
                print_vector(x32)
                print("F(x):")
                print_vector(Ft3(x32))

                k33, x33 = newton_multivar(Ft3, Jft3, [1 for i in range(80)], 1e-5, 50)
                print("\n\nn = 80")
                print("Número de iterações: %d" %k33)
                print("x:")
                print_vector(x33)
                print("F(x):")
                print_vector(Ft3(x33))

        if menu == 2:
            arquivo = input("digite o nome do arquivo: ")
            infoSat = read_file(arquivo)

            #========================================================================================================
            # Definindo F, J e funções auxiliares para calculo das coordenadas pelo GPS
            def ri(x,i):
                '''
                    vetor , indice -> numero
                    - calcula ri
                    - i deve ser a linha de r a ser referida
                '''
                tam = len(x)
                r_i = 0

                # note que var se refere a qual variavel estamos lidando: x, y, z ou w
                for var in range(tam):
                    if var < tam-1:
                        r_i += (x[var] - infoSat[i][var])**2
                    else:
                        r_i += -(infoSat[i][var] - x[var])**2
                return r_i

            def ri_linha(x,i,d):
                '''
                    vetor , indice, indice -> numero
                    - calcula ri_linha
                    - i deve ser a linha de r a ser referida
                    - d deve ser o indice da variavel que
                        estamos fazendo a derivada parcial
                '''
                tam = len(x)
                r_i_linha = 0

                if d == tam-1:
                    r_i_linha = 2*( infoSat[i][d] - x[d] )
                else:
                    r_i_linha = 2*( x[d] - infoSat[i][d] )
                
                return r_i_linha
                

            def F(x):
                tam = len(x)
                satelites = len(infoSat)
                f = [ 0 for i in range(0,tam,1)]

                # O indice d se refere a qual derivada estamos tomando.
                # O indice i se refere a qual linha estamos somando da grad g.
                for d in range(0,tam,1):
                    for i in range(0,satelites,1):
                        if d < tam-1:
                            f[d] += 4*( ri(x,i) )*( x[d] - infoSat[i][d] )
                        else:
                            f[d] += 4*( ri(x,i) )*( infoSat[i][d] - x[d] )

                return f

            def J(x):
                tam = len(x)
                satelites = len(infoSat)
                j = [[ 0 for i in range(0,tam,1)] for j in range(0,tam,1)]

                # O indice L se refere a linha da matriz jacobiana.
                # O indice c se refere a coluna da matriz jacobiana,
                #                       note que as derivadas estão conforme a linha
                #                       pois se referem a um dos termos da grad g
                # O indice i se refere ao termo do somatorio, ou seja,
                #                       linha da matriz infoSat, ou efetivamente,
                #                       a qual satelite estamos nos referindo.
                for L in range(0,tam,1):
                    for c in range(0,tam,1):
                        if c == L:
                            aux = 1
                        else:
                            aux = 0
                        for i in range(0,satelites,1):
                            if L < tam - 1:
                                j[L][c] += 4*( ri_linha(x,i,c)*(x[L] - infoSat[i][L]) + (ri(x,i)*aux) )
                            else:
                                j[L][c] += 4*( ri_linha(x,i,c)*(infoSat[i][L] - x[L]) - (ri(x,i)*aux) )
                                
                
                return j

            def F_inicial():
                v = len(infoSat[0])
                linhas = len(infoSat)
                n = linhas-1
                fini = [[ 0 for j in range(0,v,1)] for i in range(0,linhas-1,1)]
                b = [ 0 for i in range(linhas-1)]

                for i in range(0,linhas-1,1):
                    
                    for k in range(0,v,1):
                        if k < v-1:
                            fini[i][k] = 2*(infoSat[n][k] - infoSat[i][k])
                            b[i] += -infoSat[i][k]**2 + infoSat[n][k]**2

                        else:
                            fini[i][k] = 2*(infoSat[i][k] - infoSat[n][k])
                            b[i] += infoSat[i][k]**2 - infoSat[n][k]**2

                return fini,b

            #================================================================================================================

            print("Dados de entrada:")
            print_matrix(infoSat)

            A, b = F_inicial()
            At = matrixTransp(A)
            AtA = matrixmult(At,A)
            Atb = matrixMultCol(At,b)
            x = solve_LU(AtA,Atb)
            
            print("Aproximação Inicial:")
            print_vector(x)
            
            k,x = newton_multivar(F, J, x, 1e-5, 50)
            
            print("Número de iterações: %d" %k)
            print("Critério de parada: Máximo do passo em módulo dos elementos de x ser menor que 10^(-5)", end="\n\n")
            print("Coordenadas do receptor em termos de (x,y,z,w)")
            print_vector(x)
            
            c = 299792458
            T = x[3]/c

            r = vec_norm(x)
            theta = math.degrees(math.atan2(x[1],x[0]))
            phi = 90 - math.degrees(math.atan2(vec_norm([x[0],x[1]]),x[2]))
            
            print("Raio: %.3f Km" %(r/1000))
            print("Coordenadas esféricas: %.5fº N, %.5fº E" %(phi, theta))
            print("Erro de sincronização: %.3f us" %(T*1000000))

        sair = bool(int(input("Sair: 1 Continuar:  0\n")))
    
    return

# ======================================================================

def read_file(file):
    '''
        Lê arquivo
    '''
    f = open(file)
    arquivo = [x.strip('\n') for x in f.readlines()]
    arquivo = [x.strip(' ') for x in arquivo]
    linhas = int(arquivo[0])
    aux = [arquivo[i].split("  ") for i in range(1,linhas+1,1)]
    sist = []
    for i in range(0,len(aux),1):
        sist.append([])
        for j in range(0,len(aux[i]),1):
            if len(aux[i][j]) != 0:
                sist[i].append(float(aux[i][j]))
    return sist

# ======================================================================


def newton_multivar(f, jacob, x, eps, maxiter):

    k = 0
    convergiu = False
    while (not convergiu):
        
        # Verifica se a Newton falhou
        if k >= maxiter:
            return(-1,x)

        # Aplica o passo da Newton
        A = jacob(x)
        b = prod_esc(-1,f(x))
        c = solve_LU(A,b)

        # Caso o incremento do vetor x tenha norma muito pequena, consideramos que a função convergiu
        if max([abs(c[i]) for i in range(len(c))]) < eps:
            convergiu = True

        # Atualiza vetor solução
        x = vec_soma(x,c)

        # Adiciona 1 em contador de iterações
        k += 1
        
    return (k, x)
# ======================================================================
def Ft1(x):
    ft1 = [(2*(x[0]-2)),(2*(x[1]-3))]
    return ft1
def Jft1(x):
    jft1 = [[2,0],[0,2]]
    return jft1
# ======================================================================
def Ft2(x):
    ft2 = [ (4*x[0]-x[1]+x[2]-x[0]*x[3]) , (-x[0]+3*x[1]-2*x[2]-x[1]*x[3]) , (x[0]-2*x[1]+3*x[2]-x[2]*x[3]) , (x[0]**2+x[1]**2+x[2]**2-1) ]
    return ft2
def Jft2(x):
    jft2 = [ [4-x[3],-1,1,-x[0]] , [-1,3-x[3],-2,-x[1]] , [1,-2,3-x[3],-x[2]] , [2*x[0],2*x[1],2*x[3],0] ]
    return jft2
# ======================================================================
def Ft3(x):
    n = len(x) + 1
    ft3 = [ -x[i-1] +2*x[i] -x[i+1] -math.exp(x[i])/(n**2) if(i>0 and i<n-2) else 2*x[i] -x[i+1] -math.exp(x[i])/(n**2) if i==0 else -x[i-1] +2*x[i] -math.exp(x[i])/(n**2) for i in range(0,n-1,1)]
    return ft3
def Jft3(x):
    n = len(x) + 1
    jft3 = [ [ (2-math.exp(x[i])/(n**2)) if i==j else -1 if abs(j-i)==1 else 0 for j in range(n-1) ] for i in range(n-1)]
    return jft3

# ======================================================================
def matrixmult (A, B):
    rows_A = len(A)
    cols_A = len(A[0])
    rows_B = len(B)
    cols_B = len(B[0])

    if cols_A != rows_B:
      print ("Cannot multiply the two matrices. Incorrect dimensions.")
      return

    # Create the result matrix
    # Dimensions would be rows_A x cols_B
    C = [[0 for row in range(cols_B)] for col in range(rows_A)]

    for i in range(rows_A):
        for j in range(cols_B):
            for k in range(cols_A):
                C[i][j] += A[i][k] * B[k][j]
    return C

# ======================================================================

def matrixMultCol(A,b):
    '''
        Mulplica uma matriz por um vetor 
    '''
    n = len(A)
    m = len(A[0])
    Ab = [0 for i in range(0,n,1)]

    for i in range(n):
        for j in range(m):
            Ab[i] += A[i][j]*b[j]
     
    return Ab

# ======================================================================

def matrixTransp(A):
    '''
        Transpoe matriz
    '''
    m = len(A)
    n = len(A[0])

    At = [ [ A[i][j] for i in range(m) ] for j in range(n)]

    return At

# ======================================================================

def prod_esc (k, x):
    ''' numero, vetor -> vetor
        Multiplica vetor por escalar
    '''
    kx = [0 for i in range(0,len(x),1)]
    for j in range ( 0, len(x), 1 ):
        kx[j] = k*x[j]
    return kx

# ======================================================================

def vec_soma (x, y):
    ''' vetor, vetor, inteiro-> vetor
        Soma dois vetores
    '''
    t = len(x)
    soma = [0 for i in range(0,t,1)]
    for j in range ( 0, t, 1 ):
        soma[j] = x[j] + y[j]
    return soma

# ======================================================================

def vec_norm (x):
    ''' vetor -> float
        Calcula norma de um vetor
    '''
    norm_quad = 0
    for j in range ( 0, len(x), 1 ):
        norm_quad += x[j]**2
    norm = math.sqrt(norm_quad)
    return float(norm)

# ======================================================================

def print_matrix(matriz):
    ''' matriz -> null '''
    for k in range(15):
        print("=", end="")
    print("",end="\n")
    for n in range (0,len(matriz),1):
        print(" ", end="")
        for m in range (0,len(matriz[0]),1):
            print ("| %.5f" %(matriz[n][m]), end="")
        print("|", end="\n\n")
    for k in range(15):
        print("=", end="")
    print("", end="\n\n")
    return

#=======================================================================

def print_vector(vetor):
    ''' matriz -> null '''
    for k in range(15):
        print("=", end="")
    print("",end="\n")
    for n in range (0,len(vetor),1):
        print ("| %.5f" %(vetor[n]), end="")
    print("|")
    for k in range(15):
        print("=", end="")
    print("", end="\n\n")
    return

#=======================================================================

def solve_sys_U (matriz, b):
    ''' matriz, vetor -> vetor
        Resolve um sistema triangular superior por subtituição
    '''
    m = len(matriz)
    x= [0 for i in range(0,m,1)]
    for j in range (m-1,-1,-1):
        subtrair = 0
        for k in range (m-1,j,-1):
            subtrair += matriz[j][k]*x[k]
        x[j] = (b[j] - subtrair)/matriz[j][j]
        
    return x
        

# ======================================================================

def solve_sys_L (matriz, b):
    ''' matriz, vetor -> vetor
        Resolve um sistema triangular inferior por subtituição
    '''
    m = len(matriz)
    x = [0 for i in range(0,m,1)]
    for j in range (0,m,1):
        subtrair = 0
        for k in range (0,j,1):
            subtrair += matriz[j][k]*x[k]
        x[j] = (b[j] - subtrair)
        
    return x

# ======================================================================

def troca_linha(troca,p):
    '''
        Troca as linhas de um vetor "troca" de acordo com os indices do vetor "p"
    '''
    for k in range(len(troca)):
        if k != p[k]:
            guarda = troca[k]
            troca[k] = troca[p[k]]
            troca[p[k]] = guarda
    return

#=======================================================================

def solve_LU (matriz, b):
    ''' matriz, vetor -> vetor
        Resolve sistema linear por LU
    '''

    
    tam = len(matriz)
    # transpondo para nao alterar valores de matriz e b
    A = [[matriz[i][j] for j in range(tam)] for i in range(tam)]
    b_lin = [b[i] for i in range(tam)]

    # decomposicao LU
    A,p = LU_decomp(A)
    
    # troca elementos de b_lin para que fique de acordo com a decomp LU
    troca_linha(b_lin,p)
    

    # Faz Ly=b_linha
    y = solve_sys_L(A,b_lin)

    # Faz Ux=y
    x = solve_sys_U(A,y)

    
    return x
        

# ======================================================================

def LU_decomp(A):
    '''
        Faz decomposição LU de uma matriz A
    '''
    n = len(A)
    p = [0 for i in range(0,n,1)]
    for k in range(n):
        for i in range(k,n,1):
            A[i][k] = A[i][k] - sum( A[i][j]*A[j][k] for j in range(0,k,1) )
        
        # determina pivo
        l = 0
        maximo = 0
        for i in range(k,n,1):
            if maximo < abs(A[i][k]):  
                l = i
                maximo = abs(A[i][k])
        p[k] = l
        
        if k != p[k]:
            guarda = A[k]
            A[k] = A[p[k]]
            A[p[k]] = guarda
        
        for j in range(k+1,n,1):
            A[k][j] = A[k][j] - sum( A[k][i]*A[i][j] for i in range(0,k,1) )
            A[j][k] = A[j][k]/A[k][k]

    return A,p

main()
