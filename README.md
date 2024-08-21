# Implementação de Algoritmos de Grafos
## Alunos Bernardo Pavani & Anwar Hermuche

Este projeto contém a implementação de diversos algoritmos de grafos utilizando a classe `Graph` em Python para o trabalho prático da Universidade Federal de Lavras, da disciplina de Algoritmos em Grafos. A classe suporta grafos direcionados e não-direcionados, permitindo realizar operações comuns, como verificar conectividade, bipartição, existência de ciclos, e mais.

## Classe Graph

### Inicialização
```python
Graph(num_vertices, is_directed)
```
- *num_vertices*: Número de vértices no grafo.
- *is_directed*: Booleano indicando se o grafo é dirigido (True) ou não (False).

## Métodos Principais

### Adicionar Aresta: add_edge(edge_id, u, v, weight)
Adiciona uma aresta entre os vértices u e v com peso weight.

### Verificar Conectividade: is_connected()
Retorna "1" se o grafo é conexo, caso contrário "0".

### Verificar Bipartição: is_bipartite()
Retorna "1" se o grafo é bipartido, caso contrário "0".

### Verificar Ciclos: has_cycle()
Retorna "1" se o grafo possui ciclo, caso contrário "0".

### Árvore DFS: dfs_tree()
Retorna os IDs das arestas na árvore DFS a partir do vértice 0.

### Árvore BFS: bfs_tree()
Retorna os IDs das arestas na árvore BFS a partir do vértice 0.

### Ordenação Topológica: topological_sort()
Retorna a ordenação topológica dos vértices (para grafos direcionados).

### Valor do Caminho Mínimo: path_value()
Retorna o valor do caminho mais curto do vértice 0 ao vértice n-1.

### Fluxo Máximo: maximum_flow()
Retorna o valor do fluxo máximo do vértice 0 ao vértice n-1 (para grafos direcionados).

Entre outros...

## Exemplo de uso
```python
g = Graph(5, True)
g.add_edge(0, 0, 1, 10)
g.add_edge(1, 1, 2, 5)
print(g.is_connected())  # Output: "1"
print(g.topological_sort())  # Output: "0 1 2 3 4"
```

