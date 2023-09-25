import numpy as np

from qutip import *

def integrate(N, h, Jx, Jy, Jz, psi0, tlist):
    id = qeye(2)
    sz = sigmaz()
    sx = sigmax()
    sy = sigmay()

    def pauliI(N,i,pauli):
        temp = [id for x in range(N)]
        temp[i] = pauli
        return tensor(temp)

    
    H0 = 0
    for x in range(N):
        H0 += h[x] * pauliI(N,x,sz)

    H0 = H0 * (-0.5)
    H1 = 0 
    
    for x in range(N-1):
        H1 += Jx[x] * pauliI(N,x,sx) * pauliI(N,x+1,sx)
        H1 += Jy[x] * pauliI(N,x,sy) * pauliI(N,x+1,sy)
        H1 += Jz[x] * pauliI(N,x,sz) * pauliI(N,x+1,sz)
    H1 = H1 * -.5

    def H1_coeff(t, args):
        return np.exp(-(t/5.)**2)

    H = [H0, [H1, H1_coeff]]
    # evolve and calculate expectation values
    result = mesolve(H, psi0, tlist)
    return result.states


# For local test
if __name__ == "__main__":  
    N = 3            # number of spins
    h  = 1.0 * 2 * np.pi * np.ones(N) 
    Jz = 0.1 * 2 * np.pi * np.ones(N)
    Jx = 0.1 * 2 * np.pi * np.ones(N)
    Jy = 0.1 * 2 * np.pi * np.ones(N)

    # intial state, first spin in state |1>, the rest in state |0>
    psi_list = []
    psi_list.append(basis(2,1))
    for n in range(N-1):
        psi_list.append(basis(2,0))
    psi0 = tensor(psi_list)
    point = 5
    tlist = np.linspace(0, 2, point)
    sz_expt = integrate(N, h, Jx, Jy, Jz, psi0, tlist)
    print(sz_expt[point-1].full())

# For local test, your program should print

# [[ 0.        +0.j        ] 
#  [-0.49820154+0.1954353j ] 
#  [ 0.36553555+0.53850419j] 
#  [ 0.        +0.j        ] 
#  [ 0.5017984 +0.19543705j] 
#  [ 0.        +0.j        ] 
#  [ 0.        +0.j        ] 
#  [ 0.        +0.j        ]]



