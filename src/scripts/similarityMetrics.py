from functools import reduce
from nltk.metrics.segmentation import windowdiff
import segalign
import segeval

from helpers import (
    are_segments_overlapping, mass_to_sgm, massToStr, massToBinStr, segment_jaccard, seg_to_set, set_intersect_ratio
)





# This function is not needed but you may find it useful, it showcases how our metric, windowDiff, and B, behave on interesting examples
def compare_metric_behavior():
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
        a = segalign.alignment_index(h1,h2)
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
    a = segalign.alignment_index(h1,h2)
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
    a = segalign.alignment_index(h1,h2)
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
    a = segalign.alignment_index(h1,h2)
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
    a = segalign.alignment_index(h1,h2)
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
    a = segalign.alignment_index(h1,h2)
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
    a = segalign.alignment_index(h1,h2)
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
    a = segalign.alignment_index(h1,h2)
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
    a = segalign.alignment_index(h1,h2)
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

    a = segalign.alignment_index(h1,h2)
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

    a = segalign.alignment_index(h1,h2)
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
    a = segalign.alignment_index(h1,h2)
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
    a = segalign.alignment_index(h1,h2)
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
    a = segalign.alignment_index(h1,h2)
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
    a = segalign.alignment_index(h1,h2)
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
    a = segalign.alignment_index(h1,h2)
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
        a = segalign.alignment_index(h1,h2)
        b = segeval.boundary_similarity(h1,h2)
        wd = 1-windowdiff(massToBinStr(h1),massToBinStr(h2),k)
        
        print(f"z: {a:0.2f}; b: {b:0.2f}; 1-wd: {wd:0.2f};")
        print()


            

