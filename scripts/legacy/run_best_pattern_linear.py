import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "code"))

from best_pattern_search_linear import search_best

# C3 x C7, k=5
search_best(3, 7, 5)
