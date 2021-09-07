import pygame
from numpy import zeros
N       = 128
iter    = 16
t       = 0
white   = (255, 255, 255)
black   = (  0,   0,   0)

def getIndex(i,j):
    if i > N-1 : i = N-1
    elif i < 0 : i = 0
    if j > N-1 : j = N-1
    elif j < 0 : j = 0
    return i + j*N

class Fluid:
    def __init__(self, dt, diffusion, viscosity):
        self.size   = N
        self.dt     = dt
        self.diff   = diffusion
        self.visc   = viscosity

        # self.s      = zeros(N*N)
        # self.density= zeros(N*N)
        # self.Vx     = zeros(N*N)
        # self.Vy     = zeros(N*N)
        #
        # self.Vx0    = zeros(N*N)
        # self.Vy0    = zeros(N*N)

        self.s      = [0 for i in range(N*N)]
        self.density= [0 for i in range(N*N)]
        self.Vx     = [0 for i in range(N*N)]
        self.Vy     = [0 for i in range(N*N)]

        self.Vx0    = [0 for i in range(N*N)]
        self.Vy0    = [0 for i in range(N*N)]

    def fluidCubeStep(self):
        n       = self.size
        visc    = self.visc
        diff    = self.diff
        dt      = self.dt
        Vx      = self.Vx
        Vy      = self.Vy
        Vx0     = self.Vx0
        Vy0     = self.Vy0
        s       = self.s
        density = self.density

        diffuse(1, Vx0, Vx, visc, dt)
        diffuse(2, Vy0, Vy, visc, dt)

        project(Vx0, Vy0, Vx, Vy)

        advect(1, Vx, Vx0, Vx0, Vy0, dt)
        advect(2, Vy, Vy0, Vx0, Vy0, dt)

        project(Vx, Vy, Vx0, Vy0)

        diffuse(0, s, density, diff, dt)
        advect(0, density, s, Vx, Vy, dt)

    def addDensity(self, x, y, amount):
        self.density[getIndex(x, y)]    += amount

    def addVelocity(self, x, y, amountX, amountY):
        index = getIndex(x, y)
        self.Vx[index] += amountX
        self.Vy[index] += amountY

    def renderD(self, screen, scale):
        surf = pygame.Surface((scale, scale))
        surf.fill(white)
        for i in range(N):
            for j in range(N):
                x   = i*scale
                y   = j*scale
                d   = self.density[getIndex(i, j)]
                surf.set_alpha(d)
                screen.blit(surf, (x, y))

    def renderV(self, screen, scale):
        for i in range(N):
            for j in range(N):
                x   = i*scale
                y   = j*scale
                vx  = self.Vx[getIndex(i, j)]
                vy  = self.Vy[getIndex(i, j)]

                if not (abs(vx) < 0.1 and abs(vy) < 0.1):
                    pygame.draw.line(x, y, x+vx*scale, y+vy*scale)

    def fadeD(self, screen, scale):
        for d in self.density:
            if d < 0.02 : d = 0
            elif d > 254.98 : d = 255

def set_bnd(b, x):
    for i in range(1, N-1):
        if b==2:
            x[getIndex(i, 0  )]    =-x[getIndex(i, 1  )]
            x[getIndex(i, N-1)]    =-x[getIndex(i, N-2)]
        else:
            x[getIndex(i, 0  )]    = x[getIndex(i, 1  )]
            x[getIndex(i, N-1)]    = x[getIndex(i, N-2)]
    for j in range(1, N-1):
        if b==1:
            x[getIndex(0  , j)]      =-x[getIndex(1, j)]
            x[getIndex(N-1, j)]      =-x[getIndex(1, j)]
        else:
            x[getIndex(0  , j)]      = x[getIndex(1, j)]
            x[getIndex(N-1, j)]      = x[getIndex(1, j)]
    x[getIndex(0, 0)]      = 0.5 * (x[getIndex(1,   0  )] + x[getIndex(0  , 1  )])
    x[getIndex(0, N-1)]    = 0.5 * (x[getIndex(1,   N-1)] + x[getIndex(0  , N-2)])
    x[getIndex(N-1, 0)]    = 0.5 * (x[getIndex(N-2, 0  )] + x[getIndex(N-1, 1  )])
    x[getIndex(N-1, N-1)]  = 0.5 * (x[getIndex(N-2, N-1)] + x[getIndex(N-1, N-2)])



def lin_solve(b, x, x0, a, c):
    cRecip  = 1./c
    for k in range(iter):
        for j in range(1, N-1):
            for i in range(1, N-1):
                x[getIndex(i,j)]    = x0[getIndex(i,j)]
                temp                = x[getIndex(i-1, j)]
                temp               += x[getIndex(i+1, j)]
                temp               += x[getIndex(i, j+1)]
                temp               += x[getIndex(i, j-1)]
                temp               *= a
                x[getIndex(i,j)]   += temp*cRecip
    set_bnd(b, x)

def diffuse(b, x, x0, diff, dt):
    a = dt * diff * (N - 2) * (N - 2)
    lin_solve(b, x, x0, a, 1+4*a)

def project(velocX, velocY, p, div):
    for j in range(1, N-1):
        for i in range(1, N-1):
            temp    = velocX[getIndex(i+1, j)]
            temp   +=-velocX[getIndex(i-1, j)]
            temp   += velocY[getIndex(i, j+1)]
            temp   +=-velocY[getIndex(i, j-1)]
            temp   *= -0.5/N
            div[getIndex(i, j)]     = temp
            p[getIndex(i,j)]        = 0
    set_bnd(0, div)
    set_bnd(0, p)
    lin_solve(0, p, div, 1, 4)

    for j in range(1, N-1):
        for i in range(1, N-1):
            velocX[getIndex(i,j)] -= 0.5* (p[getIndex(i+1,j)] - p[getIndex(i-1,j)])*N
            velocY[getIndex(i,j)] -= 0.5* (p[getIndex(i,j+1)] - p[getIndex(i,j-1)])*N
    set_bnd(1, velocX)
    set_bnd(2, velocY)

def advect(b, d, d0, velocX, velocY, dt):
    dtx = dt * (N-2)
    dty = dtx
    Nfloat  = N
    for j in range(1, N-1):
        jfloat  = j
        for i in range(1, N-1):
            ifloat  = i
            tmp1    = dtx * velocX[getIndex(i,j)]
            tmp2    = dty * velocY[getIndex(i,j)]
            x       = ifloat - tmp1
            y       = jfloat - tmp2

            if x < 0.5 : x = 0.5
            if x > Nfloat + 0.5 : x = Nfloat + 0.5
            i0 = int(x)
            i1 = i0 + 1
            if y < 0.5 : y = 0.5
            if y > Nfloat + 0.5 : y = Nfloat + 0.5
            j0 = int(y)
            j1 = j0 + 1

            s1 = x  - i0
            s0 = 1. - s1
            t1 = y  - j0
            t0 = 1. - t1

            i0i = int(i0)
            i1i = int(i1)
            j0i = int(j0)
            j1i = int(j1)

            d[getIndex(i,j)]   = s0 * (t0 * d0[getIndex(i0i, j0i)] + t1 * d0[getIndex(i0i, j1i)])
            d[getIndex(i,j)]  += s1 * (t0 * d0[getIndex(i1i, j0i)] + t1 * d0[getIndex(i1i, j1i)])
    set_bnd(b, d)
