# Trabalho de Grafos - Professor Mayron Moreira
# Alunos: Bernardo Coelho Pavani Marinho (202310874) & Anwar Martins Hermuche (202220165)
# Turma dos dois alunos: 10A

from collections import defaultdict, deque
import heapq

class Graph:
    def __init__(self, num_vertices, is_directed):
        self.num_vertices = num_vertices
        self.is_directed = is_directed
        self.adjacency_list = defaultdict(list)
        self.edge_list = []
        self.edge_ids = {}
        
    def add_edge(self, edge_id, u, v, weight):
        self.adjacency_list[u].append((v, weight))
        self.edge_list.append((edge_id, u, v, weight))
        self.edge_ids[(u, v)] = edge_id
        if not self.is_directed:
            self.adjacency_list[v].append((u, weight))

    # método 0 - verificar se é conexo
    def is_connected(self):
      visited = set()

      # escolhendo vértice inicial
      start_vertex = 0

      if self.is_directed:
          # criando uma versão não direcionada do grafo para checar a conectividade fraca
          undirected_graph = Graph(self.num_vertices, False)
          for u in self.adjacency_list:
              for v, weight in self.adjacency_list[u]:
                  undirected_graph.add_edge(self.edge_ids[(u, v)], u, v, weight)
          
          undirected_graph._dfs_connected(start_vertex, visited)
      else:
          self._dfs_connected(start_vertex, visited)

      # conferindo se todos os vertices foram visitados
      return str(int(len(visited) == self.num_vertices))

    def _dfs_connected(self, vertex, visited):
        visited.add(vertex)
        for neighbor, _ in self.adjacency_list[vertex]:
            if neighbor not in visited:
                self._dfs_connected(neighbor, visited)

    def _reverse_graph(self):
        reversed_graph = Graph(self.num_vertices, self.is_directed)
        for u in self.adjacency_list:
            for v, weight in self.adjacency_list[u]:
                reversed_graph.add_edge(self.edge_ids[(u, v)], v, u, weight)
        return reversed_graph

    # método 1 - verificar se é bipartido
    def is_bipartite(self):
        color = [-1] * self.num_vertices  # -1 significa não colorido
        
        for start_vertex in range(self.num_vertices):
            if color[start_vertex] == -1:  # começa a BFS desse vértice
                if not self._bfs_check_bipartite(start_vertex, color):
                    return "0"  # não bipartido
        
        return "1"  # bipartido

    def _bfs_check_bipartite(self, start_vertex, color):
        queue = deque([start_vertex])
        color[start_vertex] = 0  # começando a coloração com cor 0
        
        while queue:
            u = queue.popleft()
            
            for v, _ in self.adjacency_list[u]:
                if color[v] == -1:  # se não é colorido, coloque uma cor diferente
                    color[v] = 1 - color[u]
                    queue.append(v)
                elif color[v] == color[u]:  # se os vértices adjacentes possuem a mesma cor
                    return False
        
        return True

    # método 2 - verificar se é euleriano
    def is_eulerian(self):
        if self.is_directed:
            return self._is_eulerian_directed()
        else:
            return self._is_eulerian_undirected()

    def _is_eulerian_undirected(self):
        # confere se o grafo é conexo usando o método is_connected criado anteriormente
        if self.is_connected() == "0":
            return "0"  # não Euleriano

        # vendo se todos os vértices possuem grau par 
        if all(len(self.adjacency_list[v]) % 2 == 0 for v in range(self.num_vertices)):
            return "1"  # euleriano

        return "0"  # não Euleriano

    def _is_eulerian_directed(self):
        in_degree = [0] * self.num_vertices
        out_degree = [0] * self.num_vertices

        for u in range(self.num_vertices):
            for v, _ in self.adjacency_list[u]:
                out_degree[u] += 1
                in_degree[v] += 1

        # verificando se o grau de entrada é igual ao grau de saída dos vértices (pseudossimetria)
        if all(in_degree[v] == out_degree[v] for v in range(self.num_vertices)):
            # verificando se o grafo é fortemente conexo usando o método já criado
            if self._is_strongly_connected():
                return "1"  # euleriano
            else:
                return "0"  # não Euleriano

        return "0"  # não Euleriano

    def _is_strongly_connected(self):
        # verificando a conectividade do grafo original
        if self.is_connected() == "0":
            return False
        
        # vendo se os vértices são alcançáveis no grafo reverso
        reversed_graph = self._reverse_graph()
        return reversed_graph.is_connected() == "1"

    # método 3 - verificar se possui ciclo
    def has_cycle(self):
        if self.is_directed:
            return self._has_cycle_directed()
        else:
            return self._has_cycle_undirected()

    def _has_cycle_undirected(self):
        visited = [False] * self.num_vertices
        for vertex in range(self.num_vertices):
            if not visited[vertex]:
                if self._dfs_cycle_undirected(vertex, visited, -1):
                    return "1"  # ciclo encontrado
        return "0"  # não há ciclo

    def _dfs_cycle_undirected(self, vertex, visited, parent):
        visited[vertex] = True
        for neighbor, _ in self.adjacency_list[vertex]:
            if not visited[neighbor]:
                if self._dfs_cycle_undirected(neighbor, visited, vertex):
                    return True
            elif neighbor != parent:
                return True
        return False

    def _has_cycle_directed(self):
        visited = [False] * self.num_vertices
        recursion_stack = [False] * self.num_vertices
        for vertex in range(self.num_vertices):
            if not visited[vertex]:
                if self._dfs_cycle_directed(vertex, visited, recursion_stack):
                    return "1"  # ciclo encontrado
        return "0"  # não há ciclo

    def _dfs_cycle_directed(self, vertex, visited, recursion_stack):
        visited[vertex] = True
        recursion_stack[vertex] = True

        for neighbor, _ in self.adjacency_list[vertex]:
            if not visited[neighbor]:
                if self._dfs_cycle_directed(neighbor, visited, recursion_stack):
                    return True
            elif recursion_stack[neighbor]:
                return True

        recursion_stack[vertex] = False
        return False

    # método 4 - número de componentes conexas em um grafo não orientado
    def number_of_connected_components(self):
      if self.is_directed:
        return "-1"

      visited = [False] * self.num_vertices
      component_count = 0

      for vertex in range(self.num_vertices):
          if not visited[vertex] and len(self.adjacency_list[vertex]) > 0:
            # começa uma nova DFS, significando um novo componente
            self._dfs_connected_components(vertex, visited)
            component_count += 1

      return str(component_count)

    def _dfs_connected_components(self, vertex, visited):
      visited[vertex] = True
      for neighbor, _ in self.adjacency_list[vertex]:
        if not visited[neighbor]:
          self._dfs_connected_components(neighbor, visited)

    # método 5 - número de componentes fortemente conexas em um grafo orientado
    def number_of_strongly_connected_components(self):
        if not self.is_directed:
            return "-1"

        # fazendo a primeira BFS e salvando a ordem final
        visited = [False] * self.num_vertices
        finish_stack = []

        for vertex in range(self.num_vertices):
            if not visited[vertex]:
                self._dfs_first_pass(vertex, visited, finish_stack)

        # fazendo o grafo reverso
        reversed_graph = self._reverse_graph()

        # fazendo o DFS pela segunda vez reduzindo o tempo final
        visited = [False] * self.num_vertices
        scc_count = 0

        while finish_stack:
            vertex = finish_stack.pop()
            if not visited[vertex]:
                reversed_graph._dfs_second_pass(vertex, visited)
                scc_count += 1

        return str(scc_count)

    def _dfs_first_pass(self, vertex, visited, finish_stack):
        visited[vertex] = True
        for neighbor, _ in self.adjacency_list[vertex]:
            if not visited[neighbor]:
                self._dfs_first_pass(neighbor, visited, finish_stack)
        finish_stack.append(vertex)

    def _dfs_second_pass(self, vertex, visited):
        visited[vertex] = True
        for neighbor, _ in self.adjacency_list[vertex]:
            if not visited[neighbor]:
                self._dfs_second_pass(neighbor, visited)

    # método 6 - listar os vértices de articulação
    def find_articulation_points(self):
        if self.is_directed:
            return "-1"

        visited = [False] * self.num_vertices
        disc = [-1] * self.num_vertices  # tempos de descoberta dos vértices visitados
        low = [-1] * self.num_vertices  # vértice mais baixo alcançável
        parent = [-1] * self.num_vertices  # vértices pais na árvore DFS
        articulation_points = set()  # usar um set para evitar duplicatas
        time = [0]  # inteiro mutável para rastrear o tempo da DFS

        for vertex in range(self.num_vertices):
            if not visited[vertex]:
                self._dfs_articulation(vertex, visited, disc, low, parent, articulation_points, time)

        # converter para uma lista ordenada para ordem lexicográfica e retornar como uma string separada por espaços
        articulation_points = sorted(articulation_points)
        retorno = " ".join(map(str, articulation_points))
        if retorno:
            return retorno
        else:
            return "-1"

    def _dfs_articulation(self, u, visited, disc, low, parent, articulation_points, time):
        children = 0  # número de filhos na árvore DFS
        visited[u] = True
        disc[u] = low[u] = time[0]
        time[0] += 1

        for v, _ in sorted(self.adjacency_list[u]):  # processar vizinhos em ordem lexicográfica
            if not visited[v]:
                parent[v] = u
                children += 1
                self._dfs_articulation(v, visited, disc, low, parent, articulation_points, time)

                # verificar se a subárvore enraizada em v tem uma conexão de volta para um dos ancestrais de u
                low[u] = min(low[u], low[v])

                # caso 1: u é a raiz e tem mais de um filho (considerando a condição correta para a raiz)
                if parent[u] == -1 and children > 1:
                    articulation_points.add(u)

                # caso 2: u não é raiz e o valor low de um dos seus filhos é maior ou igual ao disc[u]
                if parent[u] != -1 and low[v] >= disc[u]:
                    articulation_points.add(u)

            elif v != parent[u]:  # atualizar o valor low de u para as chamadas de função do pai
                low[u] = min(low[u], disc[v])

    # método 7 - número de arestas ponte em um grafo não direcionado
    def number_of_bridges(self):
        if self.is_directed:
          return "-1"

        visited = [False] * self.num_vertices
        disc = [-1] * self.num_vertices  # tempo de descoberta dos vértices visitados
        low = [-1] * self.num_vertices  # vértice visitado mais cedo (menor valor)
        parent = [-1] * self.num_vertices  # vértices pai na árvore DFS
        bridges_count = 0
        time = [0]  # inteiro mutável para trackear o tempo da DFS

        for vertex in range(self.num_vertices):
            if not visited[vertex]:
                bridges_count += self._dfs_bridges(vertex, visited, disc, low, parent, time)

        return str(bridges_count)

    def _dfs_bridges(self, u, visited, disc, low, parent, time):
        visited[u] = True
        disc[u] = low[u] = time[0]
        time[0] += 1
        bridges_count = 0

        for v, _ in self.adjacency_list[u]:
            if not visited[v]:
                parent[v] = u
                bridges_count += self._dfs_bridges(v, visited, disc, low, parent, time)

                # vendo se a sub árvore com raiz em v tem uma conexão para um dos ancestrais de u
                low[u] = min(low[u], low[v])

                # se o menor vértice atingível da sub árvore de v é abaixo de u na árvore DFS, então u-v é uma ponte
                if low[v] > disc[u]:
                    bridges_count += 1

            elif v != parent[u]: 
                low[u] = min(low[u], disc[v])

        return bridges_count

    # método 8 - árvore DFS
    def dfs_tree(self):
        visited = [False] * self.num_vertices
        tree_edges = []
        
        # começa no vértice 0
        self._dfs_tree(0, visited, tree_edges)
        
        # retorna a lista dos ids das arestas na árvore DFS
        return " ".join(map(str, tree_edges)) + " "

    def _dfs_tree(self, u, visited, tree_edges):
        visited[u] = True
        
        # ordena os vizinhos em ordem lexicográfica e explora cada um
        for v, _ in sorted(self.adjacency_list[u]):
            if not visited[v]:
                # Achando o id da aresta para cada aresta (u, v)
                edge_id = self.edge_ids[(u, v)] if (u, v) in self.edge_ids else self.edge_ids[(v, u)]
                tree_edges.append(edge_id)
                self._dfs_tree(v, visited, tree_edges)

    # método 9 - árvore BFS
    def bfs_tree(self):
        visited = [False] * self.num_vertices
        tree_edges = []
        queue = deque([0])  # começa o BFS no índice 0
        visited[0] = True

        while queue:
            u = queue.popleft()
            # explora os vizinhos em ordem lexicográfica
            for v, _ in sorted(self.adjacency_list[u]):
                if not visited[v]:
                    visited[v] = True
                    # encontra o id das arestas (u, v)
                    edge_id = self.edge_ids[(u, v)] if (u, v) in self.edge_ids else self.edge_ids[(v, u)]
                    tree_edges.append(edge_id)
                    queue.append(v)

        return " ".join(map(str, tree_edges)) + " "

    # método 10 - valor final da Árvore Geradora Mínima
    def mst_final_value(self):
      if self.is_directed:
          return "-1"

      # ordena as arestas pelo peso
      sorted_edges = sorted(self.edge_list, key=lambda edge: edge[3])

      # inicializa a estrutura union-find
      parent = list(range(self.num_vertices))
      rank = [0] * self.num_vertices

      def find(u):
          if parent[u] != u:
              parent[u] = find(parent[u])
          return parent[u]

      def union(u, v):
          root_u = find(u)
          root_v = find(v)
          if root_u != root_v:
              if rank[root_u] > rank[root_v]:
                  parent[root_v] = root_u
              elif rank[root_u] < rank[root_v]:
                  parent[root_u] = root_v
              else:
                  parent[root_v] = root_u
                  rank[root_u] += 1

      # aplica Kruskal
      mst_value = 0
      edges_used = 0

      for edge_id, u, v, weight in sorted_edges:
          if find(u) != find(v):
              union(u, v)
              mst_value += weight
              edges_used += 1

          # pare antes se tivermos usado n-1 arestas
          if edges_used == self.num_vertices - 1:
              break

      return str(mst_value)
    
    # método 11 - imprimir vértices em ordem topológica
    def topological_sort(self):
        if not self.is_directed:
            return "-1"

        # calcula o grau de entrada para cada vértice
        in_degree = [0] * self.num_vertices
        for u in range(self.num_vertices):
            for v, _ in self.adjacency_list[u]:
                in_degree[v] += 1

        # inicializa um min-heap com todos os vértices com grau de entrada 0
        min_heap = []
        for i in range(self.num_vertices):
            if in_degree[i] == 0:
                heapq.heappush(min_heap, i)

        topological_order = []

        # aplica o algoritmo de Kahn (com ordem lexicográfica)
        while min_heap:
            u = heapq.heappop(min_heap)
            topological_order.append(u)

            # reduz o grau de entrada em 1 unidade para cada vizinho de u
            for v, _ in self.adjacency_list[u]:
                in_degree[v] -= 1
                if in_degree[v] == 0:
                    heapq.heappush(min_heap, v)

        # a ordenação topológica inclui todos os vértices?
        if len(topological_order) != self.num_vertices:
            return "-1"  # O grafo possui ciclo

        return " ".join(map(str, topological_order))
    
    # método 12 - valor do caminho mínimo entre dois vértices
    def path_value(self):
        if not self.is_directed:
            # checando se todos os pesos são iguais
            if self._all_weights_equal():
                return "-1"

        # usando dijkstra 
        return self._dijkstra(0, self.num_vertices - 1)

    def _all_weights_equal(self):
        if not self.edge_list:
            return False
        first_weight = self.edge_list[0][3]
        for _, _, _, weight in self.edge_list:
            if weight != first_weight:
                return False
        return True

    def _dijkstra(self, start, goal):
        distances = [float('inf')] * self.num_vertices
        distances[start] = 0
        priority_queue = [(0, start)]  # (distância, vértice)

        while priority_queue:
            current_dist, u = heapq.heappop(priority_queue)

            if u == goal:
                return str(current_dist)

            if current_dist > distances[u]:
                continue

            for v, weight in self.adjacency_list[u]:
                distance = current_dist + weight

                if distance < distances[v]:
                    distances[v] = distance
                    heapq.heappush(priority_queue, (distance, v))

        return "-1"  # se nao existe caminho entre o inicio e o objetivo

    # método 13 - valor do fluxo máximo
    def maximum_flow(self):
        if not self.is_directed:
            return "-1"

        return self._edmonds_karp(0, self.num_vertices - 1)

    def _bfs_flow(self, residual_graph, source, sink, parent):
        visited = [False] * self.num_vertices
        queue = deque([source])
        visited[source] = True

        while queue:
            u = queue.popleft()

            for v in residual_graph[u]:
                if not visited[v] and residual_graph[u][v] > 0:  # capacidade residual positiva
                    queue.append(v)
                    visited[v] = True
                    parent[v] = u

                    if v == sink:
                        return True

        return False

    def _edmonds_karp(self, source, sink):
        # inicializa o grafo residual com as capacidades originais
        residual_graph = defaultdict(lambda: defaultdict(int))
        for edge_id, u, v, capacity in self.edge_list:
            residual_graph[u][v] = capacity
            residual_graph[v][u] = 0  # aresta reversa com capacidade inicial 0

        parent = [-1] * self.num_vertices
        max_flow = 0

        while self._bfs_flow(residual_graph, source, sink, parent):
            # achando o fluxo máximo através do caminho encontrado pela BFS
            path_flow = float('inf')
            s = sink

            while s != source:
                path_flow = min(path_flow, residual_graph[parent[s]][s])
                s = parent[s]

            # atualizando a capacidade residual das arestas e vértices reversos
            v = sink
            while v != source:
                u = parent[v]
                residual_graph[u][v] -= path_flow
                residual_graph[v][u] += path_flow
                v = parent[v]

            max_flow += path_flow

        return str(max_flow)

    # método 14 - fecho transitivo
    def transitive_closure(self):
        if not self.is_directed:
            return "-1"

        # inicializa o conjunto de visitados e começa a DFS no vértice 0
        reachable = set()
        visited = [False] * self.num_vertices

        # inicia a DFS no vértice 0
        self._dfs_transitive_closure(0, visited, reachable)

        # retorna os vértices alcançáveis ordenados em ordem lexicográfica
        return " ".join(map(str, sorted(reachable)))

    def _dfs_transitive_closure(self, u, visited, reachable):
        # marque o nó atual como visitado e o adiciona ao conjunto de alcançáveis 
        visited[u] = True
        reachable.add(u)

        # vizinhos em ordem lexicográfica
        for v, _ in sorted(self.adjacency_list[u]):
            if not visited[v]:
                self._dfs_transitive_closure(v, visited, reachable)

        # tentanto explorar outros vértices depois de processar todos os vértices alcançáveis de u
        for i in range(self.num_vertices):
            if not visited[i]:
                self._dfs_transitive_closure(i, visited, reachable)

    # executa os métodos das funções
    def execute_methods(self, functions_to_use):
        mapping = {0: self.is_connected, 1: self.is_bipartite, 2: self.is_eulerian, 3: self.has_cycle, 4: self.number_of_connected_components,
                      5: self.number_of_strongly_connected_components, 6: self.find_articulation_points, 7: self.number_of_bridges, 8: self.dfs_tree,
                      9: self.bfs_tree, 10: self.mst_final_value, 11: self.topological_sort, 12: self.path_value, 13: self.maximum_flow,
                      14: self.transitive_closure}
        for index in functions_to_use:
            print(mapping[index]())

def read_graph_input():
    # input
    functions_to_use = list(map(int, input().strip().split()))
    num_vertices, num_edges = map(int, input().strip().split())
    graph_type = input().strip()
    
    is_directed = (graph_type == 'direcionado')
    
    graph = Graph(num_vertices, is_directed)
    
    for _ in range(num_edges):
        edge_id, u, v, weight = map(int, input().strip().split())
        graph.add_edge(edge_id, u, v, weight)
        
    return functions_to_use, graph

functions_to_use, graph = read_graph_input()
graph.execute_methods(functions_to_use)
