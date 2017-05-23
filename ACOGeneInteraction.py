import itertools
import bisect
import random
from copy import copy

class WorldGI:
    """Nodes and edges of a particular gene network.

    :class:`WorldGI` is created from a correleation 'matrix', the nodes (genes) are represented by the column names (that 
    has to be a unique identifyer) and the length between two nodes is represented by the correlation value between
    the two of them. Additionally, :class:`WorldGI` has a variable `times` (default is 1) to determine how many times 
    the ant has to go thorough each node (gene).
    
    :class:`WorldGI` objects convert the actual nodes into node IDs, since solving does not rely on the actual 
    data in the nodes. These are accessible via the :attr:`nodes` property. Also, the name of the nodes are 
    accesible via the :attr:`nodes_names`.
    
    :class:`Edge`\\s are accessible in much the same way, except two node IDs must be passed to 
    the :func:`data` method to indicate which nodes start and end the :class:`Edge`
    """

    uid = 0
    def __init__(self, matrix, times = 1):
        self.uid = self.__class__.uid
        self.__class__.uid += 1
        self.matrix = 1 - abs(matrix)
        self.edges = self.create_edges()
        self.times = times
    
    @property
    def nodess(self):
        nodess = []
        for i in list(self.matrix.columns):
            x = y = i
            nodess.append((x, y))
        return nodess
    
    @property
    def node_names(self):
        """Node's Name"""
        nodename = []
        for i in list(self.matrix.columns):
            nodename.append(i)
        return nodename
    
    @property
    def nodes(self):
        """Node's ID"""
        return list(range(len(self.matrix.columns)))

    def create_edges(self):
        """Create edges from the nodes and the correlation matrix"""
        edges = {}
        for m in list(range(len(self.matrix.columns))):
            for n in list(range(len(self.matrix.index))):
                a, b = m, n
                edge = Edge(a, b, length = self.matrix.iloc[m,n])
                edges[m, n] = edge
        return edges

    def reset_pheromone(self, level=0.01):
        """Reset the amount of pheromone on every edge to some base *level*.
        Each time a new set of solutions is to be found, the amount of
        pheromone on every edge should be equalized to ensure un-biased initial
        conditions.(default=0.01)
        """

        for edge in self.edges.values():
            edge.pheromone = level
            
    def data(self, idx, idy=None):
        """Return the node data of a single id or the edge data of two ids."""
        try:
            if idy is None:
                return self.nodess[idx]
            else:
                return self.edges[idx, idy]
        except IndexError:
            return None


class Edge:
    """This class represents the link between starting and ending nodes.
    In addition to *start* and *end* nodes, every :class:`Edge` has *length*
    and *pheromone* properties. *length* represents the static, *a priori*
    information, whereas *pheromone* level represents the dynamic, *a
    posteriori* information.
    
    :param node start: the node at the start of the :class:`Edge`
    :param node end: the node at the end of the :class:`Edge`
    :param float length: the length of the :class:`Edge` (default=1)
    :param float pheromone: the amount of pheromone on
    the :class:`Edge`(default=0.1)
    """

    def __init__(self, start, end, length=None, pheromone=None):
        self.start = start
        self.end = end
        self.length = 1 if length is None else length
        self.pheromone = 0.1 if pheromone is None else pheromone
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        return False

class Ant:
    """
    A single independent finder of solutions to a :class:`World`.

    Each :class:`Ant` finds a solution to a world one move at a time.  They
    also represent the solution they find, and are capable of reporting which
    nodes and edges they visited, in what order they were visited, and the
    total length of the solution.

    Two properties govern the decisions each :class:`Ant` makes while finding
    a solution: *alpha* and *beta*. *alpha* controls the importance placed on
    pheromone while *beta* controls the importance placed on distance. In 
    general, *beta* should be greater than *alpha* for best results.

    :class:`Ant`\\s also have a *uid* property that can be used to identify a
    particular instance."""
    
    uid = 0
    
    def __init__(self, alpha=1, beta=3):
        """Create a new Ant for the given world."""
        self.uid = self.__class__.uid
        self.__class__.uid += 1
        self.world = None
        self.alpha = alpha
        self.beta = beta
        self.start = None
        self.distance = 0
        self.visited = []
        self.unvisited = []
        self.traveled = []

    def initialize(self, world, start=None):
        """Reset everything so that a new solution can be found."""
        self.world = world
        if start is None:
            self.start = random.randrange(len(self.world.nodes))
        else:
            self.start = start
        self.distance = 0
        self.visited = [self.start]
        self.unvisited = [n for n in self.world.nodes if n != self.start] + (self.world.times-1)*[n for n in self.world.nodes]
        self.traveled = []
        return self

    def clone(self):
        """Return a shallow copy with a new UID."""
        ant = Ant(self.alpha, self.beta)
        ant.world = self.world
        ant.start = self.start
        ant.visited = self.visited[:]
        ant.unvisited = self.unvisited[:]
        ant.traveled = self.traveled[:]
        ant.distance = self.distance
        return ant

    @property
    def node(self):
        """Most recently visited node."""
        try:
            return self.visited[-1]
        except IndexError:
            return None

    @property
    def tour(self):
        """Nodes visited by the :class:`Ant` in order."""
        return [self.world.data(i) for i in self.visited]
    
    @property
    def nodes_visited(self):
        """Nodes visited in order in the tour. """
        x = [self.world.data(i) for i in self.visited]
        y = [i for i,j in x]
        return y
    
    @property
    def path(self):
        """Edges traveled by the :class:`Ant` in order."""
        return [edge for edge in self.traveled]
        
    @property
    def tour_interactions(self):
        """ Unique interactions present in the tour."""
        c = [(x[1],y[0]) for x,y in zip(self.tour,self.tour[1:])] #edges solution tour
        return (c)
        
    @property
    def interactions(self):
        """ All the interations present in the tour."""
        c = [(x[1],y[0]) for x,y in zip(self.tour,self.tour[1:])] #edges solution tour
        h= []
        for i, j in c:
            x = i,j
            if i > j: 
                y = j,i
                h.append(y)
            if j > i:
                h.append(x)
        d = set(h)
        return (d)

    def __eq__(self, other):
        """Return ``True`` if the distance is equal to the other distance."""
        return self.distance == other.distance

    def __lt__(self, other):
        """Return ``True`` if the distance is less than the other distance."""
        return self.distance < other.distance

    def can_move(self):
        """Return ``True`` if there are moves that have not yet been made."""
        # This is only true after we have made the move back to the starting
        # node.
        return len(self.traveled) != len(self.visited)

    def move(self):
        """Choose, make, and return a move from the remaining moves."""
        remaining = [n for n in self.unvisited if n != self.node]
        choice = self.choose_move(remaining)
        return self.make_move(choice)

    def choose_move(self, choices):
        """Choose a move from all possible moves."""
        if len(choices) == 0:
            return None
        if len(choices) == 1:
            return choices[0]
        
        # Find the weight of the edges that take us to each of the choices.
        weights = []
        for move in choices:
            if move != self.node:
                edge = self.world.edges[self.node, move]
                weights.append(self.weigh(edge))
        
        # Choose one of them using a weighted probability.
        total = sum(weights)
        cumdist = list(itertools.accumulate(weights)) + [total]
        return choices[bisect.bisect(cumdist, random.random() * total)]

    def make_move(self, dest):
        """Move to the *dest* node and return the edge traveled."""
        # Since self.node simply refers to self.visited[-1], which will be
        # changed before we return to calling code, store a reference now.
        ori = self.node
        # When dest is None, all nodes have been visited but we may not
        # have returned to the node on which we started. If we have, then
        # just do nothing and return None. Otherwise, set the dest to the
        # node on which we started and don't try to move it from unvisited
        # to visited because it was the first one to be moved.
        if dest is None:
            if self.can_move() is False:
                return None
            dest = self.start   # last move is back to the start
        else:
            self.visited.append(dest)
            self.unvisited.remove(dest)
        
        edge = self.world.edges[ori, dest]
        self.traveled.append(edge)
        self.distance += edge.length
        return edge

    def weigh(self, edge):
        """Calculate the weight of the given *edge*."""
        pre = 1 / (edge.length or 1)
        post = edge.pheromone
        return post ** self.alpha * pre ** self.beta


class SolverGI:
    """This class contains the functionality for finding one or more solutions
    for a given :class:`WorldGI`.
    
    :param float alpha: relative importance of pheromone (default=1)
    :param float beta: relative importance of distance (default=3)
    :param float rho: percent evaporation of pheromone (0..1, default=0.8)
    :param float q: total pheromone deposited by each :class:`Ant` after
                    each iteration is complete (>0, default=1)
    :param float t0: initial pheromone level along each :class:`Edge` of the
                     :class:`WorldGI` (>0, default=0.01)
    :param int limit: number of iterations to perform (default=100)
    :param float ant_count: how many :class:`Ant`\\s will be used 
                            (default=10)
    :param float elite: multiplier of the pheromone deposited by the elite
                        :class:`Ant` (default=0.5)
    """

    def __init__(self, alpha = 1, beta = 3, rho = 0.8, q = 1, t0 = 0.01, limit = 100 , ant_count = 10, elite = 0.05):
        self.alpha = alpha
        self.beta = beta
        self.rho = rho
        self.q = q
        self.t0 = t0
        self.limit = limit
        self.ant_count = ant_count
        self.elite = elite
    
    def create_colony(self, world):
        """Create a set of :class:`Ant`\\s and initialize them to the given *world*"""
        if self.ant_count < 1:
            return self.round_robin_ants(world, len(world.nodes))
        return self.random_ants(world, self.ant_count)
        
    def reset_colony(self, colony):
        """Reset the *colony* of :class:`Ant`\\s such that each :class:`Ant` is ready 
        to find a new solution"""
        for ant in colony:
            ant.initialize(ant.world)
        
    def aco(self, colony):
        """Return the best solution by performing the ACO meta-heuristic."""
        self.find_solutions(colony)
        self.global_update(colony)
        return sorted(colony)[0]
        
    def solve(self, world):
        """Return the single shortest path found through the given *world*."""
        world.reset_pheromone(self.t0)
        global_best = None
        colony = self.create_colony(world)
        for i in range(self.limit):
            self.reset_colony(colony)
            local_best = self.aco(colony)
            if global_best is None or local_best < global_best:
                global_best = copy(local_best)
            self.trace_elite(global_best)
        return global_best
    
    def solutions(self, world):
        """Return successively shorter paths through the given *world*."""
        world.reset_pheromone(self.t0)
        global_best = None
        colony = self.create_colony(world)
        for i in range(self.limit):
            self.reset_colony(colony)
            local_best = self.aco(colony)
            if global_best is None or local_best < global_best:
                global_best = copy(local_best)
                yield global_best
            self.trace_elite(global_best)
    
    def round_robin_ants(self, world, count):
        """Returns a list of :class:`Ant`\\s distributed to the nodes of the 
        world in a round-robin fashion."""
        starts = world.nodes
        n = len(starts)
        return [
            Ant(self.alpha, self.beta).initialize(
                world, start=starts[i % n])
            for i in range(count)
        ]
        
    def random_ants(self, world, count, even=False):
        """Returns a list of :class:`Ant`\\s distributed to the nodes of the 
        world in a random fashion."""
        ants = []
        starts = world.nodes
        n = len(starts)
        if even:
            # Since the caller wants an even distribution, use a round-robin 
            # method until the number of ants left to create is less than the
            # number of nodes.
            if count > n:
                for i in range(self.ant_count // n):
                    ants.extend([
                        Ant(self.alpha,self.beta).initialize(
                            world, start=starts[j])
                        for j in range(n)
                    ])
            # Now (without choosing the same node twice) choose the reamining
            # starts randomly.
            ants.extend([
                Ant(self.alpha, self.beta).initialize(
                    world, start=starts.pop(random.randrange(n - i)))
                for i in range(count % n)
            ])
        else:
            # Just pick random nodes.
            ants.extend([
                Ant(self.alpha, self.beta).initialize(
                    world, start=starts[random.randrange(n)]) 
                for i in range(count)
            ])
        return ants
    
    def find_solutions(self, ants):
        """Let each :class:`Ant` find a solution.
        Makes each :class:`Ant` move until each can no longer move."""
        # This loop occurs exactly as many times as there are ants times nodes,
        # but that is only because every ant must visit every node. It may be
        # more efficient to convert it to a counting loop... but what 
        # flexibility would we loose in terms of extending the solver class?
        ants_done = 0
        while ants_done < len(ants):
            ants_done = 0
            for ant in ants:
                if ant.can_move():
                    edge = ant.move()
                    self.local_update(edge)
                else:
                    ants_done += 1

    def local_update(self, edge):
        """Evaporate some of the pheromone on the given *edge*."""
        edge.pheromone = max(self.t0, edge.pheromone * self.rho)

    def global_update(self, ants):
        """Update the amount of pheromone on each edge according to the fitness
        of solutions that use it."""
        ants = sorted(ants)[:len(ants) // 2]
        for a in ants:
            p = self.q / a.distance
            for edge in a.path:
                edge.pheromone = max(
                    self.t0,
                    (1 - self.rho) * edge.pheromone + p)

    def trace_elite(self, ant):
        """Deposit pheromone along the path of a particular ant.
        This method is used to deposit the pheromone of the elite :class:`Ant`
        at the end of each iteration."""
        if self.elite:
            p = self.elite * self.q / ant.distance
            for edge in ant.path:
                edge.pheromone += p

