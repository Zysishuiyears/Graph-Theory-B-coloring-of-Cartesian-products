# -*- coding: utf-8 -*-
"""
Created on Sun Oct 12 22:28:17 2025

@author: JZX
"""

# import numpy as np

# for i in range(3, 13, 2):
#     length_list = []
#     for j in range(0, 20, 2):
#         length = 0
#         length = i+5*j
        
#         length_list.append(length)
        
#     print("基础大小:", i, "实际大小:", length_list)
    
from sympy import primerange

def uncovered_odd_unordered():
    primes = list(primerange(2, 101))
    a_pool = {m for p in primes for m in range(p, 101, p)} & set(range(3, 101))
    b_pool = {m for p in primes for m in range(p, 101, p)} & set(range(3, 101))
    #b_pool = set(range(5, 101, 5)) 
    
    covered = {frozenset([a, b]) for a in a_pool for b in b_pool if a % 2 and b % 2}
    all_odds = {frozenset([x, y]) for x in range(3, 101, 2) for y in range(3, 101, 2)}
    return all_odds - covered, len(all_odds)


if __name__ == "__main__":
    uncovered, total = uncovered_odd_unordered()
    # 排序：先转 tuple 再 sorted
    ordered = [tuple(sorted(fs)) for fs in uncovered]
    ordered.sort(key=lambda t: (t[0], t[1] if len(t) == 2 else t[0]))

    # 写文件
    with open("uncovered_odd_ordered.txt", "w") as f:
        for t in ordered:
            if len(t) == 1:                      # 自对 (x,x)
                f.write(f"{t[0]},{t[0]}\n")
            else:                                # (x,y) 且 x≤y
                f.write(f"{t[0]},{t[1]}\n")

    print(f"总无序奇数对：{total}，未被覆盖：{len(uncovered)}，已写入 uncovered_odd_ordered.txt")