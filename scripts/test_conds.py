#!/usr/bin/env python3
"""
plot_conductivities.py

Leitura simples de um arquivo .msh (formato Gmsh ASCII com blocos $Nodes e $NodeData)
e plot 3D dos valores de condutividade apenas para nós cuja condutividade != 0.3805.

Uso:
    python plot_conductivities.py /caminho/para/conductivities_00001.msh

Dependências:
    numpy, matplotlib

Instalação:
    pip install numpy matplotlib
"""

import sys
import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 import necessário para 3D

def parse_msh_nodes_and_nodedata(path):
    """
    Retorna (nodes_dict, nodedata_dict)
    nodes_dict: {node_id: (x,y,z)}
    nodedata_dict: {node_id: value}
    O parser é tolerante a espaços em branco e pula tags dentro de $NodeData.
    """
    nodes = {}
    nodedata = {}

    with open(path, 'r', encoding='utf-8', errors='ignore') as f:
        lines = [ln.rstrip('\n') for ln in f]

    i = 0
    L = len(lines)
    while i < L:
        line = lines[i].strip()
        if line == "$Nodes":
            # Next line: number of nodes
            i += 1
            if i >= L:
                break
            try:
                n_nodes = int(lines[i].strip())
            except ValueError:
                # fallback: if nodes block uses a different format, try to skip
                n_nodes = 0
            i += 1
            for _ in range(n_nodes):
                if i >= L:
                    break
                parts = lines[i].split()
                if len(parts) >= 4:
                    try:
                        nid = int(parts[0])
                        x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                        nodes[nid] = (x, y, z)
                    except ValueError:
                        # ignore malformed
                        pass
                i += 1
            # skip until $EndNodes if any mismatch
            while i < L and lines[i].strip() != "$EndNodes":
                i += 1
            # move past $EndNodes
        elif line == "$NodeData":
            # Gmsh $NodeData format:
            # next line: number-of-string-tags (ns)
            # next ns lines: strings
            # next line: number-of-real-tags (nr)
            # next nr lines: reals
            # next line: number-of-int-tags (ni)
            # next ni lines: ints (often one of them is the number of entries)
            # After tags: lines with "node-number value" for each node data entry
            i += 1
            if i >= L:
                break
            # parse string tags
            try:
                ns = int(lines[i].strip())
            except ValueError:
                ns = 0
            i += 1
            for _ in range(ns):
                if i < L:
                    # skip string tag lines
                    i += 1

            # parse real tags
            if i < L:
                try:
                    nr = int(lines[i].strip())
                except ValueError:
                    nr = 0
                i += 1
            else:
                nr = 0
            for _ in range(nr):
                if i < L:
                    i += 1

            # parse int tags
            if i < L:
                try:
                    ni = int(lines[i].strip())
                except ValueError:
                    ni = 0
                i += 1
            else:
                ni = 0
            for _ in range(ni):
                if i < L:
                    # some int tags may include the number of entries
                    # but we'll not rely on that strictly; we skip them
                    i += 1

            # now read data lines until $EndNodeData
            while i < L and lines[i].strip() != "$EndNodeData":
                ln = lines[i].strip()
                if ln == "":
                    i += 1
                    continue
                # each data line is usually: node_id value [possibly more values]
                toks = ln.split()
                if len(toks) >= 2:
                    try:
                        nid = int(toks[0])
                        val = float(toks[1])
                        nodedata[nid] = val
                    except ValueError:
                        # sometimes tags or non-data lines appear; try to skip
                        pass
                i += 1
            # move past $EndNodeData (it will be skipped by outer loop increment)
        i += 1

    return nodes, nodedata

def plot_nodes(nodes, nodedata, ignore_value=0.3805, tol=1e-12, show=True, save_path=None):
    """
    Plota os nós cujo valor em nodedata difere de ignore_value (com tolerância tol).
    """
    xs = []
    ys = []
    zs = []
    vals = []

    for nid, coord in nodes.items():
        if nid in nodedata:
            val = nodedata[nid]
            if not math.isclose(val, ignore_value, rel_tol=0.0, abs_tol=tol):
                x, y, z = coord
                xs.append(x)
                ys.append(y)
                zs.append(z)
                vals.append(val)

    if len(xs) == 0:
        print("Nenhum nó com condutividade diferente de", ignore_value)
        return

    xs = np.array(xs)
    ys = np.array(ys)
    zs = np.array(zs)
    vals = np.array(vals)

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')
    sc = ax.scatter(xs, ys, zs, c=vals, cmap='viridis', s=12, depthshade=True)
    cbar = fig.colorbar(sc, ax=ax, shrink=0.6)
    cbar.set_label('Condutividade')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f'Nós com condutividade != {ignore_value} (count={len(xs)})')

    # Ajuste de aspecto aproximado
    try:
        max_range = np.array([xs.max()-xs.min(), ys.max()-ys.min(), zs.max()-zs.min()]).max() / 2.0
        mid_x = (xs.max()+xs.min()) * 0.5
        mid_y = (ys.max()+ys.min()) * 0.5
        mid_z = (zs.max()+zs.min()) * 0.5
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)
    except Exception:
        pass

    if save_path:
        fig.savefig(save_path, dpi=200)
        print(f"Figura salva em: {save_path}")

    if show:
        plt.show()
    else:
        plt.close(fig)

def main():
    if len(sys.argv) < 2:
        print("Uso: python plot_conductivities.py /caminho/para/conductivities_00001.msh")
        sys.exit(1)

    path = Path(sys.argv[1])
    if not path.exists():
        print("Arquivo não encontrado:", path)
        sys.exit(1)

    print("Lendo arquivo:", path)
    nodes, nodedata = parse_msh_nodes_and_nodedata(path)
    print(f"Nós lidos: {len(nodes)}; entradas em NodeData: {len(nodedata)}")
    import pdb;pdb.set_trace()
    # Ajuste do valor a ignorar (padrão 0.3805). Pode ser modificado aqui.
    ignore_val = 0.3805
    # Plot (mostra janela interativa). Se quiser salvar automaticamente, passe save_path.
    plot_nodes(nodes, nodedata, ignore_value=ignore_val, tol=1e-12, show=True, save_path=None)

if __name__ == "__main__":
    main()