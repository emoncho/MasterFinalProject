import pandas
import ACOGeneInteractions

Example =pd.DataFrame([(1,    0.7647, 0.1982, 0.8013,0.8018,0.3838,0.1912, 0.4326),
                (0.7647, 1,     0.5101, 0.9538,0.9603,0.2135,0.6497,0.6267),
                (0.1982, 0.5101, 1,     0.5962,0.4584,0.2996, 0.4551,0.4270),
                (0.8013, 0.9538, 0.5962, 1,     0.9779,0.4535, 0.5668,0.6465),
                (0.8018, 0.9603, 0.4584, 0.9779, 1  ,  0.4009, 0.5796,0.5855),
                (0.3838, 0.2135, 0.2996, 0.4535,0.4009,1    ,-0.0175,0.3807),
                (0.1912, 0.6497, 0.4551, 0.5668,0.5796,-0.0175,1    ,0.5159),
                (0.4326, 0.6267, 0.4270, 0.6465,0.5855,0.3807,0.5159,1  )
                ],  columns=['uvrD','lexA','umuDC','recA','uvrA','uvrY','ruvA','polB'],
                index=['uvrD','lexA','umuDC','recA','uvrA','uvrY','ruvA','polB'])

world = WorldGI(Example, times = 1)
world.nodes
world.node_names

solver1 = SolverGI()
solution = solver.solve(world)
print(solution.tour_interactions) #all the interactions
print(solution.interactions)      #unique interactions

world = WorldGI(Example, times = 2)
world.nodes
world.node_names
solver = SolverGI()
solution = solver.solve(world)
print(solution.tour_interactions) #all the interactions
print(solution.interactions)      #unique interactions

world = WorldGI(Example, times = 3)
world.nodes
world.node_names

solver = SolverGI()
solution = solver.solve(world)
print(solution.tour_interactions) #all the interactions
print(solution.interactions)      #unique interactions