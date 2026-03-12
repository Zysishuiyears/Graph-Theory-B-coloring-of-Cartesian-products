import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "code"))

from best_pattern_search import search_best

m, n, k = 3, 7, 5
(best_HV, info) = search_best(m, n, k, max_solutions=10000)

print("\n=== 最优结果 ===")
print("Score   =", info[0])
print("Periods =", info[1])
print("Penalty =", info[2])
print("方案文本 =", info[3])
