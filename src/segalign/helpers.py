from typing import Tuple, List
#This file contains multiple helper functions to convert between segmentation formats and calculate frequently needed info


def set_intersect_ratio(a,b):
    """Returns the directed intersect ratio from a to b"""
    return len(a & b) / len(a)


#returns the length of a segment of the form [start,end]
def segLen(seg1: Tuple[int,int]):
    [i,j] = seg1
    return j-i + 1


def seg_to_set(seg1: Tuple[int,int]):
    """Convert a segment :math:`s_i`, represented by its start and stop indices :math:`(s_j, s_k)`,
    to a Python ``set`` of the integers :math:`{j, j+1, j+2, ..., k}`"""
    return set(range(seg1[0], seg1[1] + 1))


def massToBinStr(massList):
    str = ''
    for segSize in massList:
        if segSize<1:
            raise ValueError(f"Error in mass list. All segments should have size >1. massList is {massList}")
        str += '0'*(segSize-1)+'1'*(1 if segSize>=1 else 0)
    return str[:-1]


def segDist(seg1: Tuple[int,int] ,seg2: Tuple[int,int],normFactor=1):
    [i,j] = seg1
    [k,l] = seg2

    return (abs(k-i) + abs (l-j), (abs(k-i) + abs (l-j))/normFactor)


def are_segments_overlapping(seg1: Tuple[int,int] ,seg2: Tuple[int,int]):
    """Returns True if two segments `seg1` and `seg2` overlap.

    Both segments are represented as tuples (start, stop)
    """
    if seg1 is None or seg2 is None:
        return False
    return not (seg1[1] < seg2[0] or seg2[1] < seg1[0])


def isSoftTransp(gSgmChunk: Tuple[Tuple[int,int],Tuple[int,int]], hSgmChunk: Tuple[Tuple[int,int],Tuple[int,int]]):
    g1 = gSgmChunk[0]
    g2 = gSgmChunk[1]

    h1 = hSgmChunk[0]
    h2 = hSgmChunk[1]

    j1 = segmentJaccard(g1,h1)
    j2 = segmentJaccard(g2,h2)

    return j1 and j2

def segment_jaccard(seg1: Tuple[int,int] ,seg2: Tuple[int,int]):
    """Compute jaccard index for two segments, ``seg1`` and ``seg2``"""
    i, j  = seg1
    k, l  = seg2

    #no overlap
    if j < k or l < i:
        return 0

    # r1 is subset of r2
    if k <= i and j <= l:
        union_size = l - k + 1
        intersect_size = union_size - (i - k) - (l - j)

    # r2 is subset of r1
    if i <= k and l <= j:
        union_size = j - i + 1
        intersect_size = union_size - (k - i) - (j - l)
    
    # r1, then r2 with intersect
    if i < k and k <= j and j < l:
        union_size = l - i + 1
        intersect_size = j - k + 1

    # r2, then r1 with intersect
    if k<i and i<=l and l<j:
        union_size = j-k+1
        intersect_size = l-i+1

    return intersect_size / union_size


def setJaccard(set1, set2):

    return (len(set1&set2))/(len(set1|set2)) 


def mass_to_sgm(mass_list: List[int]):
    """Convert a list of segment masses to a list of (start, stop) pairs"""
    total = 0
    seg = []
    for size in mass_list:
        seg.append((total, total + size - 1))
        total += size
    return seg


def massToStr(massList: List[int]):
    
    chars = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
    
    lastChar = len(chars)-1
    charIndex = 0
    
    sectionStrings = []
    for segSize in massList:
        
        section = []

        for _ in range(segSize):
            
            nextChar = chars[charIndex]
            section.append(nextChar)

            charIndex = 0 if charIndex == lastChar else charIndex+1
        
        sectionStrings.append('.'.join(section))
    
    return '|'.join(sectionStrings)

def writeLinesToFile(lines,fileAddress):
    with open(fileAddress,'w', encoding='utf8') as o:
        for l in lines:
            o.write(l+'\n')

def genSgms(nElems,nSegs):
    if nSegs<1 or nElems<1 or nSegs>nElems:
        raise ValueError(f"Invalid parameters. nElems is {nElems} and nSegs is {nSegs}")

    if nSegs == 1:
        return [[nElems]]
    
    r = []
    maxSegSize = nElems - (nSegs-1)
    for size in range(1,maxSegSize+1):
        segment = [size]
        for partition in genSgms(nElems-size,nSegs-1):
            segmentation = segment+partition
            r.append(segmentation)
    
    return r

def genSgmsDP(nElems,nSegs):
    if nSegs<1 or nElems<1 or nSegs>nElems:
        raise ValueError(f"Invalid parameters. nElems is {nElems} and nSegs is {nSegs}")

    #column for nSegs = 1
    baseColumn = [None] + [ [[n]] for n in range(1,nElems+1)]

    if nSegs == 1:
        return baseColumn[-1]

    prevColumn = baseColumn
    nextColumn = None
    
    for k in range(2,nSegs+1):
        nextColumn = [[None]]
        for n in range(1,nElems+1):
            r = []
            maxSegSize = n - (k - 1)
            for size in range(1,maxSegSize+1):
                segment = [size]
                for partition in prevColumn[n-size]:
                    segmentation = segment + partition
                    r.append(segmentation)
            
            nextColumn.append(r)
        
        prevColumn = nextColumn
    
    return nextColumn[-1]

def sgmGenerator(nElems):

    #column for nSegs = 1
    baseColumn = [None] + [ [[n]] for n in range(1,nElems+1)]

    prevColumn = baseColumn
    nextColumn = None
    
    for k in range(2,nElems):
        nextColumn = [[None]]
        for n in range(1,nElems+1):
            r = []
            maxSegSize = n - (k - 1)
            for size in range(1,maxSegSize+1):
                segment = [size]
                for partition in prevColumn[n-size]:
                    segmentation = segment + partition
                    r.append(segmentation)
            
            nextColumn.append(r)
        
        prevColumn = nextColumn
        yield nextColumn[-1]

def sgmPairGenerator(nElems):

    #column for nSegs = 1
    baseColumn = [None] + [ [[n]] for n in range(1,nElems+1)]

    prevColumn = baseColumn
    nextColumn = None
    
    k = 2
    nextColumn = [[None]]
    for n in range(1,nElems+1):
        r = []
        maxSegSize = n - (k - 1)
        for size in range(1,maxSegSize+1):
            segment = [size]
            for partition in prevColumn[n-size]:
                segmentation = segment + partition
                r.append(segmentation)
        
        nextColumn.append(r)
    
    prevColumn = nextColumn
    return nextColumn[-1]



