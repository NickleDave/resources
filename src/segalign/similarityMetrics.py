from math import inf
from functools import reduce
from typing import List

from nltk.metrics.segmentation import windowdiff
import segeval

from helpers import (
    are_segments_overlapping, mass_to_sgm, massToStr, massToBinStr, segment_jaccard, seg_to_set, set_intersect_ratio
)


# This is our alignment-based segmentation similarity metric; it runs in O(m1+m2).
# It receives 2 segmentations of the format [x0,x2,...x_n] where each element is a positive integer that represents the size of a segment
def alignment_index(s1: List[int], s2: List[int]):
    # ---- verify segmentations can be compared
    # get total number of elements in segmentations
    l1 = sum(s2)
    l2 = sum(s1)

    if l1 != l2:
        raise ValueError(
            f"Segmentations have different number of elements (the sum of segment lengths)."
            f"Number of elements in ``s1``: {l1}; Number of elements in ``s2``: {l2}"
        )

    # convert segmentations of the form [2,2] into the form [(0,1),(2,3)];
    # each segment is represented by the range of elems inside it (inclusive)
    s1 = mass_to_sgm(s1)
    s2 = mass_to_sgm(s2)

    # get the number of segments in each segmentation
    n = len(s1)
    m = len(s2)

    # ---- set up
    # list to store directed edges from alignment
    directed_edges = []
    # pointers to navigate s1(i) and s2(j)
    i = j = 0
    # variables to keep track of the "best" (highest intersect ratio) candidate segment to align to,
    # for both the ith segment in s1 and the jth segment in s2
    top_max_intersect = bottom_max_intersect = -inf
    top_best_bottom_candidate = bottom_best_top_candidate = None

    # ---- iterate left to right, starting at 0,0, until all segments in s1 and s2 have been aligned
    while i < n and j < m:

        # if the current ith and jth segments overlap, they are mutual candidates for alignment
        if are_segments_overlapping(s1[i], s2[j]):
            top_set = seg_to_set(s1[i])
            bottom_set = seg_to_set(s2[j])

            # get the intersect ratios aligning j to i and i to j
            top_int_ratio = set_intersect_ratio(top_set, bottom_set)
            bottom_int_ratio = set_intersect_ratio(bottom_set, top_set)

            # ---- potentially update best candidates for i and j;
            # candidates are updated only if
            # a) the intersect ratio is higher,
            # b) the intersect ratio is tied and the corresponding jaccard indexes are higher
            if (
                    top_int_ratio > top_max_intersect
            ) or (
                    top_int_ratio == top_max_intersect and 
                    segment_jaccard(s1[i], s2[j]) > segment_jaccard(s1[i], s2[top_best_bottom_candidate])
            ):
                top_max_intersect = top_int_ratio
                top_best_bottom_candidate = j

            if (
                    bottom_int_ratio > bottom_max_intersect
            ) or (
                    bottom_int_ratio == bottom_max_intersect and
                    segment_jaccard(s1[i], s2[j]) > segment_jaccard(s1[bottom_best_top_candidate], s2[j])):
                bottom_max_intersect = bottom_int_ratio
                bottom_best_top_candidate = i

        # if the ith segment in s1 reaches further right than the jth segment in s2,
        # then s2 has seen all candidates for alignment.
        # get the jaccard index with the best candidate for j so far and save the weighted edge
        # the jth segment in s2 has been aligned, so increase the j pointer and reset corresponding vars
        if s1[i][1] > s2[j][1]:
            score = segment_jaccard(s1[bottom_best_top_candidate],s2[j])
            edge = (f'b-{j}',f't-{bottom_best_top_candidate}', score)
            directed_edges.append(edge)

            j += 1
            bottom_best_top_candidate = None
            bottom_max_intersect = -inf

        # if the jth segment in s2 reaches further right than the ith segment in s1,
        # then s1 has seen all candidates for alignment.
        # get the jaccard index with the best candidate for i so far and save the weighted edge
        # the ith segment in s1 has been aligned, so increase the i pointer and reset corresponding vars
        elif s1[i][1] < s2[j][1]:
            score = segment_jaccard(s1[i], s2[top_best_bottom_candidate])
            edge = (f't-{i}',f'b-{top_best_bottom_candidate}', score)
            directed_edges.append(edge)

            i += 1
            top_best_bottom_candidate = None
            top_max_intersect = -inf

        # if the ith and jth segment end on the same element, then both have seen all candidates for alignment
        # get the jaccard index with the best candidate for i and j so far and save the weighted edges
        # i and j have been aligned, so increase pointers and reset corresponding vars
        elif s1[i][1] == s2[j][1]:
            score = segment_jaccard(s1[bottom_best_top_candidate], s2[j])
            edge = (f'b-{j}',f't-{bottom_best_top_candidate}',score)
            directed_edges.append(edge)

            j += 1
            bottom_best_top_candidate = None
            bottom_max_intersect = -inf

            score = segment_jaccard(s1[i],s2[top_best_bottom_candidate])
            edge = (f't-{i}',f'b-{top_best_bottom_candidate}',score)
            directed_edges.append(edge)

            i += 1
            top_best_bottom_candidate = None
            top_max_intersect = -inf
    
    # verify that each segment in s1 and s2 has been aligned
    assert len(directed_edges) == n + m

    # keep only undirected edges (deletes duplicates)
    undirected_edges = set((frozenset([x[0], x[1]]), x[2]) for x in directed_edges)

    #sum up and average edge weights
    sum_edge_weights = reduce(lambda x,y: x + y[1], undirected_edges, 0)
    avg = sum_edge_weights / len(undirected_edges)

    return avg


# This function is not needed but you may find it useful, it showcases how our metric, windowDiff, and B, behave on interesting examples
def compareMetricBehavior():
    print("Slide")

    h1=[1,8,1]
    h1s = massToStr(h1)
    alternates = [ [9,1], [2,7,1], [3,6,1], [4,5,1], [5,4,1], [6,3,1], [7,2,1], [8,1,1]]
    k = max(1,round(reduce(lambda x,y: x+y, h1)/len(h1)/2))
    print(f"k: {k}")
    for h2 in alternates:
        h2s = massToStr(h2)
        print(h1s)
        print(h2s)
        a = alignmentIndex(h1,h2)
        b = segeval.boundary_similarity(h1,h2)
        wd = 1-windowdiff(massToBinStr(h1),massToBinStr(h2),k)
        print(f"z: {a:0.2f}; b: {b:0.2f}; 1-wd: {wd:0.2f};")

    print("Misc-------------------------------------")

    h1=[3,3,3,3]
    h2 = [3,2,4,3]
    h1s = massToStr(h1)
    h2s = massToStr(h2)
    k = max(1,round(reduce(lambda x,y: x+y, h1)/len(h1)/2))
    print(f"k: {k}")
    print(h1s)
    print(h2s)
    a = alignmentIndex(h1,h2)
    b = segeval.boundary_similarity(h1,h2)
    wd = 1-windowdiff(massToBinStr(h1),massToBinStr(h2),k)
    print(f"z: {a:0.2f}; b: {b:0.2f}; 1-wd: {wd:0.2f};")

    h1=[3,3,3,3]
    h2 = [3,3,1,2,3]
    h1s = massToStr(h1)
    h2s = massToStr(h2)
    k = max(1,round(reduce(lambda x,y: x+y, h1)/len(h1)/2))
    print(f"k: {k}")
    print(h1s)
    print(h2s)
    a = alignmentIndex(h1,h2)
    b = segeval.boundary_similarity(h1,h2)
    wd = 1-windowdiff(massToBinStr(h1),massToBinStr(h2),k)
    print(f"z: {a:0.2f}; b: {b:0.2f}; 1-wd: {wd:0.2f};")

    h1=[3,3,3,3]
    h2 = [6,3,3]
    h1s = massToStr(h1)
    h2s = massToStr(h2)
    k = max(1,round(reduce(lambda x,y: x+y, h1)/len(h1)/2))
    print(f"k: {k}")
    print(h1s)
    print(h2s)
    a = alignmentIndex(h1,h2)
    b = segeval.boundary_similarity(h1,h2)
    wd = 1-windowdiff(massToBinStr(h1),massToBinStr(h2),k)
    print(f"z: {a:0.2f}; b: {b:0.2f}; 1-wd: {wd:0.2f};")

    h1=[3,3,3,3]
    h2 = [3,2,1,1,2,3]
    h1s = massToStr(h1)
    h2s = massToStr(h2)
    k = max(1,round(reduce(lambda x,y: x+y, h1)/len(h1)/2))
    print(f"k: {k}")
    print(h1s)
    print(h2s)
    a = alignmentIndex(h1,h2)
    b = segeval.boundary_similarity(h1,h2)
    wd = 1-windowdiff(massToBinStr(h1),massToBinStr(h2),k)
    print(f"z: {a:0.2f}; b: {b:0.2f}; 1-wd: {wd:0.2f};")



    print("Cross-Boundary Transpositions-------------------------------------")

    h1=[5,5,1,3]
    h2 = [7,3,1,3]
    h1s = massToStr(h1)
    h2s = massToStr(h2)
    k = max(1,round(reduce(lambda x,y: x+y, h1)/len(h1)/2))
    print(f"k: {k}")
    print(h1s)
    print(h2s)
    a = alignmentIndex(h1,h2)
    b = segeval.boundary_similarity(h1,h2)
    wd = 1-windowdiff(massToBinStr(h1),massToBinStr(h2),k)
    print(f"z: {a:0.2f}; b: {b:0.2f}; 1-wd: {wd:0.2f};")

    h1=[5,5,1,3]
    h2 = [5,4,1,4]
    h1s = massToStr(h1)
    h2s = massToStr(h2)
    k = max(1,round(reduce(lambda x,y: x+y, h1)/len(h1)/2))
    print(f"k: {k}")
    print(h1s)
    print(h2s)
    a = alignmentIndex(h1,h2)
    b = segeval.boundary_similarity(h1,h2)
    wd = 1-windowdiff(massToBinStr(h1),massToBinStr(h2),k)
    print(f"z: {a:0.2f}; b: {b:0.2f}; 1-wd: {wd:0.2f};")




    print("Constant Cost Transp-------------------------------------")

    h1=[2,2,5,5]
    h2 = [2,2,4,6]
    h1s = massToStr(h1)
    h2s = massToStr(h2)
    k = max(1,round(reduce(lambda x,y: x+y, h1)/len(h1)/2))
    print(f"k: {k}")
    print(h1s)
    print(h2s)
    a = alignmentIndex(h1,h2)
    b = segeval.boundary_similarity(h1,h2)
    wd = 1-windowdiff(massToBinStr(h1),massToBinStr(h2),k)
    print(f"z: {a:0.2f}; b: {b:0.2f}; 1-wd: {wd:0.2f};")

    h1=[2,2,5,5]
    h2 = [3,1,5,5]
    h1s = massToStr(h1)
    h2s = massToStr(h2)
    k = max(1,round(reduce(lambda x,y: x+y, h1)/len(h1)/2))
    print(f"k: {k}")
    print(h1s)
    print(h2s)
    a = alignmentIndex(h1,h2)
    b = segeval.boundary_similarity(h1,h2)
    wd = 1-windowdiff(massToBinStr(h1),massToBinStr(h2),k)
    print(f"z: {a:0.2f}; b: {b:0.2f}; 1-wd: {wd:0.2f};")


    print("Vanishing Transp-------------------------------------")
    h1=[8,2,1]
    h2 = [2,8,1]
    h1s = massToStr(h1)
    h2s = massToStr(h2)
    k = max(1,round(reduce(lambda x,y: x+y, h1)/len(h1)/2))
    print(f"k: {k}")
    print(h1s)
    print(h2s)

    a = alignmentIndex(h1,h2)
    b = segeval.boundary_similarity(h1,h2)
    wd = 1-windowdiff(massToBinStr(h1),massToBinStr(h2),k)
    print(f"z: {a:0.2f}; b: {b:0.2f}; 1-wd: {wd:0.2f};")

    h1=[7,3,1]
    h2 = [3,7,1]
    h1s = massToStr(h1)
    h2s = massToStr(h2)
    k = max(1,round(reduce(lambda x,y: x+y, h1)/len(h1)/2))
    print(f"k: {k}")
    print(h1s)
    print(h2s)

    a = alignmentIndex(h1,h2)
    b = segeval.boundary_similarity(h1,h2)
    wd = 1-windowdiff(massToBinStr(h1),massToBinStr(h2),k)
    print(f"z: {a:0.2f}; b: {b:0.2f}; 1-wd: {wd:0.2f};")


    h1=[9,5,1]
    h2 = [5,9,1]
    h1s = massToStr(h1)
    h2s = massToStr(h2)
    k = max(1,round(reduce(lambda x,y: x+y, h1)/len(h1)/2))
    print(f"k: {k}")
    print(h1s)
    print(h2s)
    a = alignmentIndex(h1,h2)
    b = segeval.boundary_similarity(h1,h2)
    wd = 1-windowdiff(massToBinStr(h1),massToBinStr(h2),k)
    print(f"z: {a:0.2f}; b: {b:0.2f}; 1-wd: {wd:0.2f};")

    h1=[8,6,1]
    h2 = [5,9,1]
    h1s = massToStr(h1)
    h2s = massToStr(h2)
    k = max(1,round(reduce(lambda x,y: x+y, h1)/len(h1)/2))
    print(f"k: {k}")
    print(h1s)
    print(h2s)
    a = alignmentIndex(h1,h2)
    b = segeval.boundary_similarity(h1,h2)
    wd = 1-windowdiff(massToBinStr(h1),massToBinStr(h2),k)
    print(f"z: {a:0.2f}; b: {b:0.2f}; 1-wd: {wd:0.2f};")

    h1=[8,6,1]
    h2 = [6,8,1]
    h1s = massToStr(h1)
    h2s = massToStr(h2)
    k = max(1,round(reduce(lambda x,y: x+y, h1)/len(h1)/2))
    print(f"k: {k}")
    print(h1s)
    print(h2s)
    a = alignmentIndex(h1,h2)
    b = segeval.boundary_similarity(h1,h2)
    wd = 1-windowdiff(massToBinStr(h1),massToBinStr(h2),k)
    print(f"z: {a:0.2f}; b: {b:0.2f}; 1-wd: {wd:0.2f};")

    h1=[7,7,1]
    h2 = [6,8,1]
    h1s = massToStr(h1)
    h2s = massToStr(h2)
    k = max(1,round(reduce(lambda x,y: x+y, h1)/len(h1)/2))
    print(f"k: {k}")
    print(h1s)
    print(h2s)
    a = alignmentIndex(h1,h2)
    b = segeval.boundary_similarity(h1,h2)
    wd = 1-windowdiff(massToBinStr(h1),massToBinStr(h2),k)
    print(f"z: {a:0.2f}; b: {b:0.2f}; 1-wd: {wd:0.2f};")

    h1=[7,7,1]
    h2 = [7,7,1]
    h1s = massToStr(h1)
    h2s = massToStr(h2)
    k = max(1,round(reduce(lambda x,y: x+y, h1)/len(h1)/2))
    print(f"k: {k}")
    print(h1s)
    print(h2s)
    a = alignmentIndex(h1,h2)
    b = segeval.boundary_similarity(h1,h2)
    wd = 1-windowdiff(massToBinStr(h1),massToBinStr(h2),k)
    print(f"z: {a:0.2f}; b: {b:0.2f}; 1-wd: {wd:0.2f};")


    
    print("Maximal Seg Test-------------------------------------")
    h1 = [1 for x in range(10)]
    for x in range(10):
        sum = 1*x
        h2 = [1 for y in range(x)] + [10-sum]
        s1 = massToStr(h1)
        s2 = massToStr(h2)
        k = max(1,round(reduce(lambda x,y: x+y, h1)/len(h1)/2))
        print(f"k: {k}")
        print(s1)
        print(s2)
        a = alignmentIndex(h1,h2)
        b = segeval.boundary_similarity(h1,h2)
        wd = 1-windowdiff(massToBinStr(h1),massToBinStr(h2),k)
        
        print(f"z: {a:0.2f}; b: {b:0.2f}; 1-wd: {wd:0.2f};")
        print()


            

