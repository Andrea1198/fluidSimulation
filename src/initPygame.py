import pygame
from src.functions import Fluid, N
from numpy import pi, cos, sin
from random import randint
scale   = 16
width   = scale*N
height  = scale*N
screen  = pygame.display.set_mode((width, height))

fluid = Fluid(0.02, 0, 1e-7)
def start():
    t       = 0
    pygame.init()
    pygame.display.set_caption("Fluids")
    running = True
    while running:
        screen.fill((0,0,0))
        cx  = int(0.5*width/scale)
        cy  = int(0.5*height/scale)
        for i in range(-1, 2):
            for j in range(-1, 2):
                fluid.addDensity(cx+i, cy+j, randint(50,150))
        for i in range(2):
            angle   = randint(0,100)/100*2*pi
            v       = [cos(angle)*0.02, sin(angle)*0.02]
            t      +=0.01
            fluid.addVelocity(cx, cy, v[0], v[1])
        fluid.fluidCubeStep()
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False
        fluid.renderD(screen, scale)
        pygame.display.update()
